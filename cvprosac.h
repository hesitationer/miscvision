/*M/////////////////////////////////////////////////////////////////////////////////////////
// cvRansac -- Random Sample Consensus algorithm
//
// Adapted From:
// ransac.m
// % Peter Kovesi
// % School of Computer Science & Software Engineering
// % The University of Western Australia
// % pk at csse uwa edu au    
// % http://www.csse.uwa.edu.au/~pk
// %
// % References:
// %    M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
// %    for model fitting with applications to image analysis and automated
// %    cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981
// %
// %    Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
// %    Computer Vision". pp 101-113. Cambridge University Press, 2001
/////////////////////////////////////////////////////////////////////////////////////////M*/

#include <cv.h>
#include <highgui.h>

/*//////////////////////////////////////////////////////////////////////////////////////////
// Function that solves for model parameters given data items, and a list of indices 
// @param data 		MxN matrix where the rows are data items
// @param indices   nfitx1 CV_32SC1 matrix where the rows are indicies into rows of data.  
//                  The parameter nfit is specified as an argument to cvRansac
// @param model     application specific matrix containing the output model parameters
// @return          0 if fit was successful, 1 if the fit was degenerate
/////////////////////////////////////////////////////////////////////////////////////////M*/
typedef int (* CvRansacFitFunc) (CvMat * data, CvMat * indices, CvMat * model);

/*//////////////////////////////////////////////////////////////////////////////////////////
// Function that computes distance between actual data point and model, and labels inliers
// @param data 		MxN matrix where the rows are data items
// @param model     application specific matrix containing the model parameters
// @param dist 		inlierdist passed in from cvRansac, threshold to mark inlying points
// @param inliers   application specific matrix indicating inlying data items.  This is not
//                  used internally and can be NULL 
// @return          number of data points that fit the model (i.e. number of inliers)
/////////////////////////////////////////////////////////////////////////////////////////M*/
typedef int (* CvRansacDistFunc) (CvMat * data, CvMat * model, double dist, CvMat * inliers);

/*//////////////////////////////////////////////////////////////////////////////////////////
// Function that determines if the given data points will produce a degenerate model.  Only
// use this if it is more efficient to determine if given data points are deficient than
// actually solve for the model
// @param data 		MxN matrix where the rows are data items
// @param indices   Mx1 CV_32SC1 matrix where the rows are indicies into rows of data
// @return 			0 if not degenerate, 1 if degenerate
/////////////////////////////////////////////////////////////////////////////////////////M*/
typedef int (* CvRansacDegenFunc) (CvMat * data, CvMat * indices);

/*M////////////////////////////////////////////////////////////////////////////////////////
// Performs RANSAC optimization
//  @param data       MxN matrix where the rows are data items 
//  @param inliers    matrix indicating inliers -- format, size depends on application
//  @param model      preallocated array of model parameters to be solved for -- format, size depends on application
//  @param fitfunc    function that solves for model parameters given an array of data items
//  @param distfunc   function that determines how well a data item fits a model -- must return number of inliers
//  @param degenfunc  Function that returns 1 if a set of datapoints produces a degenerate model.
//  				  0, otherwise.   This is used to discard random samples that do not
//                    result in useful models.  (Can be NULL)
//  @param maxtrials  The maximum number of times to resample
//  @param nfit       The minimum number of samples from data required by
//                    fitfunc to fit a model
//  @param inlierdist The distance threshold between a data point and the model
//                    used to decide whether the point is an inlier or not.
//  @return           Number of inliers for the model 
////////////////////////////////////////////////////////////////////////////////////////M*/
int cvProsac(CvArr * data, CvArr * dataW, CvArr * inliers, CvArr * model, CvRansacFitFunc  fitfunc, 
			 CvRansacDistFunc distfunc, CvRansacDegenFunc  degenfunc, 
			 int maxtrials, int nfit, double inlierdist, int maxdegentrials, unsigned int seed=0);
