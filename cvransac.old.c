#include "cvransac.h"
#include <cv.h>
#include <highgui.h>
#include <stdio.h>

typedef int (* CvRansacFitFunc) (CvMat * data, CvMat * indices, CvMat * model);
typedef int (* CvRansacDistFunc) (CvMat * data, CvMat * model, double dist, CvMat * inliers);
typedef int (* CvRansacDegenFunc) (CvMat * data, CvMat * indices);

int icvRansacMat(CvMat * data, CvMat * inliers, CvMat * model, CvRansacFitFunc fitfunc, 
			 CvRansacDistFunc distfunc, CvRansacDegenFunc degenfunc, 
			 int maxtrials, int nfit, double inlierdist, int maxdegentrials){
	CvMat * best_model = cvCloneMat(model);
	CvMat * tmp;
	CvMat * ind = cvCreateMat(nfit, 1, CV_32SC1);
	CvMat * best_inliers = inliers? cvCloneMat(inliers) : NULL;
	int bestscore = 0;
	
	// random number generator -- this should definitely be an argument
	CvRNG rng = cvRNG();	
	
	for(int i=0; i<maxtrials; i++){
		// Select at random s datapoints to form a trial model, M.
		// In selecting these points we have to check that they are not in
		// a degenerate configuration.
		int degenerate = 1;
		for(int j=0; j<maxdegentrials && degenerate; j++){
			// Generate random list of indices 
			cvRandArr(&rng, ind, CV_RAND_UNI, cvScalarAll(0), cvScalarAll(data->rows));
						            
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
		// of elements in x that are inliers.  Additionally, if M is a cell
		// array of possible models 'distfn' will return the model that has
		// the most inliers.  After this call M will be a non-cell object
		// representing only one model.
		int ninliers = distfunc(data, model, inlierdist, inliers);
		
		//distfunc(data, model, inlierdist, inliers, ninliers, err);

		if(ninliers > bestscore){    // Largest set of inliers so far...
			// || ninliers==bestinliers && err < besterr
			// besterr = err
			bestscore = ninliers;    // Record data for this model
			if(inliers) CV_SWAP(best_inliers, inliers, tmp);
			CV_SWAP(best_model, model, tmp);

			// Update estimate of N, the number of trials to ensure we pick, 
			// with probability p, a data set with no outliers.
			// double fracinliers =  ninliers/npts;
			// double pNoOutliers = 1 -  fracinliers^s;
			// pNoOutliers = max(eps, pNoOutliers);  // Avoid division by -Inf
			// pNoOutliers = min(1-eps, pNoOutliers); // Avoid division by 0.
			
			// Number of trials
			// N = log(1-p)/log(pNoOutliers);
		}	

		//fprintf(stdout, "trial %d out of %d         \n", i, maxtrials);
	}
	
	// copy back to output structures
	if(inliers) CV_SWAP(best_inliers, inliers, tmp);
	CV_SWAP(best_model, model, tmp);
	return bestscore;
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
			 int maxtrials, int nfit, double inlierdist, int maxdegentrials){
	// convert CvArr to CvMat
	CvMat dataM, inliersM, modelM;
	CvMat *dataMp, *inliersMp, *modelMp;
	dataMp    = cvGetMat(data, &dataM);
	inliersMp = inliers ? cvGetMat(inliers, &inliersM) : NULL;
	modelMp   = cvGetMat(model, &modelM);
	return icvRansacMat(dataMp, inliersMp, modelMp, fitfunc, distfunc, degenfunc, 
			           maxtrials, nfit, inlierdist, maxdegentrials);
}

