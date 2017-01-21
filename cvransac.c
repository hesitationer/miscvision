#include "cvransac.h"
#include <cv.h>
#include <highgui.h>
#include <cvext.h>
#include <stdio.h>
#include <limits.h>  // INT_MAX

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

/// @param  w probability of data point being 'good'
/// @param  n number of data points to form model
/// @param  p desired probability of sampling a model without outliers 
/// @return k number of trials required to ensure with probability p 
///           that n data points that contain no outliers have been selected
int icvRansacNumTrials(float w, int n, float p){
	if(w<=0) return INT_MAX;  // sort of undefined
	if(w>=1) return 0;
	// pFail = probability of selecting a model containing an outlier
	//       = 1 - w^n
	// z = probability of seeing only bad models in k trials
	// z = (pFail)^k
	// log(z) = k log(pFail)
	// z = 1 - p
	// k = log(1-p)/log(pFail)
	return log(1-p)/log(1-pow( w, n));
}

int icvRansacMat(CvMat * data, CvMat * inliers, CvMat * model, 
		CvRansacFitFunc fitfunc, CvRansacDistFunc distfunc, CvRansacDegenFunc degenfunc, 
		int maxtrials, int nfit, double inlierdist, int maxdegentrials,
		unsigned int seed = 0){

	CvMat * best_model = cvCloneMat(model);
	CvMat * tmp;
	CvMat * ind = cvCreateMat(data->rows, 1, CV_32SC1);
	CvMat * best_inliers = inliers? cvCloneMat(inliers) : NULL;
	int bestinliers = 0;
	double besterr = 0;
	int ninliers;
	double err;
	double p = .99;
	int data_size = data->rows;

	CvRNG rng = cvRNG( seed );

	// intitialize indices
	for(int i=0; i<data_size; i++){
		ind->data.i[i] = i;
	}

	for(int i=0; i<maxtrials; i++){
		// Select at random datapoints to form a trial model, M.
		// In selecting these points we have to check that they are not in
		// a degenerate configuration.
		int degenerate = 1;
		for(int j=0; j<maxdegentrials && degenerate; j++){
			// Generate random list of indices 
			// cvRandArr(&rng, ind, CV_RAND_UNI, cvScalarAll(0), cvScalarAll(data->rows));
			cvRandShuffle( ind, &rng );
						            
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

		if(ninliers > bestinliers){    // Largest/best set of inliers so far...
		  // || ninliers==bestinliers && err < besterr){
			
			besterr = err;
			bestinliers = ninliers;    
			if(inliers) CV_SWAP(best_inliers, inliers, tmp);
			CV_SWAP(best_model, model, tmp);

			// Update estimate of N, the number of trials
			double pGood =  (double)ninliers/data_size;
			maxtrials = MIN(icvRansacNumTrials(pGood, nfit, p), maxtrials);

			printf("ninliers: %d\n", ninliers);
			//printf("fracinliers = %f pNoOutliers = %f\n", fracinliers, pNoOutliers);
			printf("maxtrials updated to %d\n",  maxtrials);
		}	

		//fprintf(stdout, "trial %d out of %d\n", i, maxtrials);
	}
	
	// copy back to output structures if pointers have been switched 
	if(inliers && best_inliers != inliers) cvCopy(best_inliers, inliers);
	if(model != best_model) cvCopy(best_model, model);
	return bestinliers;
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
int cvRansac(CvArr * data, CvArr * inliers, CvArr * model, CvRansacFitFunc  fitfunc, 
			 CvRansacDistFunc distfunc, CvRansacDegenFunc  degenfunc, 
			 int maxtrials, int nfit, double inlierdist, int maxdegentrials, unsigned int seed)
{
	// convert CvArr to CvMat
	CvMat dataM, inliersM, modelM;
	CvMat *dataMp, *inliersMp, *modelMp;
	dataMp    = cvGetMat(data, &dataM);
	inliersMp = inliers ? cvGetMat(inliers, &inliersM) : NULL;
	modelMp   = cvGetMat(model, &modelM);
	return icvRansacMat(dataMp, inliersMp, modelMp, fitfunc, distfunc, degenfunc, 
			           maxtrials, nfit, inlierdist, maxdegentrials, seed);
}
