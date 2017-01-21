#include "cvransac.h"
#include <cv.h>
#include <cxmisc.h>
#include <cvext.h>  // need cvCeil1
#include <highgui.h>
#include <stdio.h>

/*typedef int (* CvRansacFitFunc) (void * data, int * indices, int nidx, void * outModel);
typedef int (* CvRansacDistFunc) (void * data, void * model, double inlierdist, 
		                          double * outErr, void * outInliers, int * outNumInliers);
typedef int (* CvRansacDegenFunc) (void * data, int * indices, int nidx);
*/

typedef int (* CvRansacFitFunc) (CvMat * data, CvMat * indices, CvMat * model);
typedef int (* CvRansacDistFunc) (CvMat * data, CvMat * model, double dist, CvMat * inliers);
typedef int (* CvRansacDegenFunc) (CvMat * data, CvMat * indices);

/*M////////////////////////////////////////////////////////////////////////////////////////
//  @model data       MxN matrix where the rows are data items 
//  @model model      preallocated array of model parameters
//  @model fitfunc    function that solves for model parameters given an array of data items
//  @model distfunc   function that determines how well a data item fits a model
//  @model degenfunc  Function that determines whether a
//                    set of datapoints will produce a degenerate model.
//                    This is used to discard random samples that do not
//                    result in useful models.  (Can be NULL)
//  @model maxtrials  The maximum number of times to resample
//  @model nfit       The minimum number of samples from data required by
//                    fitfunc to fit a model
//  @model inlierdist The distance threshold between a data point and the model
//                    used to decide whether the point is an inlier or not.
////////////////////////////////////////////////////////////////////////////////////////M*/

/*struct CvFitParam {
	void * data;
	int data_size;
	void * model;
	void * inliers;
	int nfit;
	int inlierdist;

	CvFitFunc fitfunc;
	CvDistFunc distfunc;
	CvDegenFunc degenfunc;
};

int cvFitRansac( CvFitParam * param, float prob, unsigned int seed ){ 

}
*/

#define cmp_aux( a, b ) \
	(aux[a] > aux[b])

CV_IMPLEMENT_QSORT_EX( icvSortIndices, int, cmp_aux, float *);

// T_n in the paper
// sample_size = m
// num_points = N
// max_samples = T_N
float avg_samples(int n, int sample_size, int num_points, int max_samples){
	float tn=max_samples;
	for(int i=0; i<sample_size; i++){
		tn *= (n-i);
		tn /= (num_points-i);
	}
	return tn;
}

void icvSample(int * arr, int len, int max, CvRNG * rng){
	int pop = max;
	for(int samp=len; samp>0; samp--){
		float cumprob = 1.0;
		float x = cvRandReal( rng );
		while (x<cumprob){
			cumprob -= cumprob * samp / pop;
			pop -= 1;
			assert(pop>=0);
		}
		arr[ samp-1 ] = max-pop-1;
		//printf("icvSample %d = %d\n", samp-1, arr[samp-1]);
	}
}

void icvShuffle(float * arr, int len, CvRNG * rng){
	float tmp;
	for(int i=1; i<len; i++){
		int idx = cvRandInt( rng )%(len-i)+i;
		// swap arr[i-1] and arr[idx]
		CV_SWAP(arr[i-1], arr[idx], tmp);
	}
}

float cvBinomial1(int n, int r){
	float x=1.0;
	//int k = MAX(r, n-r);
	int j = MIN(r, n-r);
	int m = MIN(n-j, j);
	for(int i=0; i<m; i++){
		x*=(n-i)/(i+1);
	}
	if(n-j<j){
		for(int i=m; i<j; i++){
			x*=1/i;
		}
	}
	else{
		for(int i=m; i<(n-j); i++){
			x*=i;
		}
	}
	return x;
}
float icvProsacPRandomSupport(float beta, int i, int m, int n){
	float p = pow(beta, i-m)*pow(1-beta, n-i+m)*cvBinomial1(n-m, i-m);
	//printf("icvProsacPRandomSupport(%f, %d, %d, %d)=%f\n", beta, i, m, n, p);
	return p;
}
int icvProsacTerminationLength(float beta, int m, int n, float psi){
	float sum=0;
	for(int i=n; i>=m; i--){
		sum+=icvProsacPRandomSupport(beta, i, m, n);
		//printf("icvProsacTerminationLength sum=%f\n", sum);
		if(sum > psi){
			return MIN(i+1,n);
		}
	}
	return m;
}

float icvProsacPUncontaminated(int ninliers, int m, int n){
	float prod = 1.0;
	for(int i=0; i<m; i++){
		prod *= (float)(ninliers - i)/(n-i);
	}
	return prod;
}

int icvProsacNumTrials( int ninliers, int m, int n, float thresh){
	if(ninliers<m) return INT_MAX;
	return cvRound( log(thresh)/log(1-icvProsacPUncontaminated(ninliers, m, n)) );
}

int icvProsacMat(CvMat * data, CvMat * dataW, CvMat * inliers, CvMat * model, 
		CvRansacFitFunc fitfunc, CvRansacDistFunc distfunc, CvRansacDegenFunc degenfunc, 
		int maxtrials, int sample_size, double inlierdist, int maxdegentrials,
		unsigned int seed = 0){

	CvMat * sorted_idx = cvCreateMat(data->rows, 1, CV_32SC1);
	CvMat * best_model = cvCloneMat(model);
	CvMat * tmp;
	CvMat * ind = cvCreateMat(sample_size, 1, CV_32SC1);
	CvMat * best_inliers = inliers? cvCloneMat(inliers) : NULL;
	float pOutlierSupported = .4;
	float pOutlierSupportedThresh = .01;
	float pMissThresh = 0.01;

	int best_ninliers = 0;
	double besterr = 0;
	int ninliers;
	double err;
	//double p = .01;

	CvRNG rng = cvRNG( seed );
	
	// PROSAC variables
	int t = 0;
	int n = sample_size;
	int num_samples_n = 1;      // T_n'
	int data_size = data->rows; // N
	int termination_length = data_size; // n*
	// sample_size =  m

	// sort indices by weights associated with data
	for(int i=0; i<sorted_idx->rows; i++){
		sorted_idx->data.i[i] = i;
	}
	icvSortIndices( sorted_idx->data.i, sorted_idx->rows, dataW->data.fl );
	//printf("data indices sorted by weight\n");
	//for(int i=0; i<sorted_idx->rows; i++){
	//	printf("%d = %f\n", sorted_idx->data.i[i], dataW->data.fl[ sorted_idx->data.i[i] ]);
	//}


	for(int i=0; i<maxtrials && best_ninliers<termination_length; i++){

		// Select at random datapoints to form a trial model, M.
		// In selecting these points we have to check that they are not in
		// a degenerate configuration.
		int degenerate = 1;
		for(int j=0; j<maxdegentrials && degenerate; j++){
			//printf("Trial %d of %d iteration %d\n", i, maxtrials, j);
			// select hypothesis generation set
			t = t+1;

			//printf("t = %d num_samples_n=%d fl(num_samples_n)=%f n=%d\n", t, num_samples_n, avg_samples(n,  sample_size, data_size, maxtrials), n);

			// t==T_n' && n<n*
			if( t==num_samples_n && n<data_size){
				// T_n+1' = T_n' + ceil( T_n+1 - T_n )
				num_samples_n += cvCeil1( avg_samples(n+1, sample_size, data_size, maxtrials) - 
						                   avg_samples(n, sample_size, data_size, maxtrials) );
				n++;
			}

			// select from the top n points
			if(num_samples_n < t){
				icvSample(ind->data.i, sample_size-1, n-1, &rng);
				ind->data.i[sample_size-1] = n-1;
			}
			else{
				icvSample(ind->data.i, sample_size, n, &rng);
			}

			// convert ranks to data indices
			for(int k=0; k<sample_size; k++){
				//printf("%d -> %d = %d\n", k, ind->data.i[k],  sorted_idx->data.i[ ind->data.i[k] ]);
				ind->data.i[k] = sorted_idx->data.i[ ind->data.i[k] ];
			}
						            
			// Test that these points are not a degenerate configuration.
			if(degenfunc){
				degenerate = degenfunc(data, ind);
			}
			else{
				degenerate=0;
			}

			if (!degenerate){

				// Fit model to this random selection of data points.
				// Note that M may represent a set of models that fit the data in
				// this case M will be a cell array of models
				degenerate = fitfunc(data, ind, model);
		
				// Depending on your problem it might be that the only way you
				// can determine whether a data set is degenerate or not is to
				// try to fit a model and see if it succeeds. 
			} 
		}
		if(degenerate){
			return 0;
		}

		// Once we are out here we should have some kind of model...        
		// Evaluate distances between points and model returning the indices
		// of elements in x that are inliers.  
		ninliers = distfunc(data, model, inlierdist, inliers);

		if(ninliers > best_ninliers){    // Largest/best set of inliers so far...
		  // || ninliers==best_ninliers && err < besterr){
			
			besterr = err;
			best_ninliers = ninliers;    
			if(inliers) CV_SWAP(best_inliers, inliers, tmp);
			CV_SWAP(best_model, model, tmp);

			termination_length = MIN( termination_length,
					icvProsacTerminationLength(pOutlierSupported, sample_size, termination_length, pOutlierSupportedThresh));
			maxtrials = MIN( icvProsacNumTrials(best_ninliers, sample_size, termination_length, pMissThresh), maxtrials );
			printf("ninliers: %d\n", ninliers);
			printf("set termination_length to %d\n", termination_length);
			printf("set maxtrials to %d\n", maxtrials);

		}	
	}
	
	// copy back to output structures if pointers have been switched 
	if(inliers && best_inliers != inliers) cvCopy(best_inliers, inliers);
	if(model != best_model) cvCopy(best_model, model);
	
	return best_ninliers;
}

/*M////////////////////////////////////////////////////////////////////////////////////////
//  @model data       MxN matrix where the rows are data items 
//  @model model      preallocated array of model parameters
//  @model fitfunc    function that solves for model parameters given an array of data items
//  @model distfunc   function that determines how well a data item fits a model
//  @model degenfunc  Function that determines whether a
//                    set of datapoints will produce a degenerate model.
//                    This is used to discard random samples that do not
//                    result in useful models.  (Can be NULL)
//  @model maxtrials  The maximum number of times to resample
//  @model nfit       The minimum number of samples from data required by
//                    fitfunc to fit a model
//  @model inlierdist The distance threshold between a data point and the model
//                    used to decide whether the point is an inlier or not.
////////////////////////////////////////////////////////////////////////////////////////M*/
int cvProsac(CvArr * data, CvArr * dataW, CvArr * inliers, CvArr * model, CvRansacFitFunc  fitfunc, 
			 CvRansacDistFunc distfunc, CvRansacDegenFunc  degenfunc, 
			 int maxtrials, int nfit, double inlierdist, int maxdegentrials, unsigned int seed)
{
	// convert CvArr to CvMat
	CvMat dataM, dataWM, inliersM, modelM;
	CvMat *dataMp, *dataWMp, *inliersMp, *modelMp;
	dataMp    = cvGetMat(data, &dataM);
	dataWMp   = cvGetMat(dataW, &dataWM);
	inliersMp = inliers ? cvGetMat(inliers, &inliersM) : NULL;
	modelMp   = cvGetMat(model, &modelM);
	return icvProsacMat(dataMp, dataWMp, inliersMp, modelMp, fitfunc, distfunc, degenfunc, 
			           maxtrials, nfit, inlierdist, maxdegentrials, seed);
}
