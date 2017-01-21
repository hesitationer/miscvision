#include <stdio.h>
#include <cv.h>
#include <highgui.h>
#include <matrix.hh>
#include "cvprosac.h"
#include <cv/circle_detector.hpp>


#if 0 
int fit_func(CvMat * data, CvMat * ind, CvMat * model){
	// fit a circle to the data points data[ind]
	// only use the first two 'white' data points
	// data[ind[0]], data[ind[1]]
	// c = midpoint between data[0], data[1]
	// r = .5 distance between data[0], data[1]
	
	// Find first two white data points:
	int * ind_data = ind->data.i;
	int i, idx1=-1, idx2=-1;
	for(i=0; i<ind->rows; i++){
		idx1 = ind_data[i]; 	
		if( cvGetReal2D(data, idx1, 2) != 0 ) break;
	}
	for( i++ ; i<ind->rows; i++){
		idx2 = ind_data[i]; 	
		if( idx2!=idx1 && cvGetReal2D(data, idx2, 2) != 0 ) break;
	}
	if(idx1==-1 || idx2==-1){
		return 1;
	}
	
	//fprintf(stdout, "indx=%d, %d\n", idx1, idx2);
	double cx,cy,r;
	double x1,x2,y1,y2;
	x1 = cvGetReal2D(data, idx1, 0);
	y1 = cvGetReal2D(data, idx1, 1);
	x2 = cvGetReal2D(data, idx2, 0);
	y2 = cvGetReal2D(data, idx2, 1);
	cx = x1 + .5*(x2-x1);
	cy = y1 + .5*(y2-y1);
	r = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) )/2;

	cvSetReal1D( model, 0, cx);
	cvSetReal1D( model, 1, cy);
	cvSetReal1D( model, 2, r);
	
	//fprintf(stdout, "circle model:
	return 0;
}
int dist_func(CvMat * data, CvMat * model, double dist, CvMat * inliers){
	double cx,cy,r;
	int ninliers=0;
	cx = cvGetReal1D(model,0);
	cy = cvGetReal1D(model,1);
	r = cvGetReal1D(model,2);
	for(int i=0; i<data->rows; i++){
		double x = cvGetReal2D(data, i, 0);
		double y = cvGetReal2D(data, i, 1);
		double c = cvGetReal2D(data, i, 2);
		int incircle = ( (x-cx)*(x-cx)+(y-cy)*(y-cy) <= r );
		if( incircle && c!=0 ||  // white pixel in circle
			!incircle && c==0 ){ // black pixel outside circle
			cvSetReal1D(inliers, i, 0);
		}
		else{
			cvSetReal1D(inliers, i, 1);
			ninliers++;
		}
	}
	return ninliers; 
}
#endif

CvCircle32f cvFitCircleRobust1( CvArr * array, CvArr * weights, CvArr * inliers, int niterations=25, double mindist=2.0){
	CvSize size = cvGetSize( array );
	Matrix<float> model(3,1);

	cvProsac( array, weights, inliers, &model, icvFitCircleRobustFitFunc, icvFitCircleRobustDistFunc, NULL, niterations, 3, mindist, 10);
	return *(CvCircle32f *) model[0];
}

int main(){
	cvRedirectError(cvSegReport);

	Matrix<float> data(100, 2);
	Matrix<float> weights(100, 1);
	Matrix<int> inliers(100,1);
	int radius = 50;
	while(1){
		for(int i=0; i<data.rows; i++){
			float theta = drand48()*M_PI*2;
			float offset = 0;
			if(drand48()>.77){
				offset=drand48()*radius;
				weights[i][0] = 1.0 - offset/radius;
			}
			else{
				weights[i][0] = drand48();
			}
			data[i][0] = cos(theta)*radius + 2*radius + offset;
			data[i][1] = sin(theta)*radius + 2*radius + offset;
		}
		IplImage * im = cvCreateImage(cvSize(500,500), 8, 3);
		cvNamedWindow("win",0);
		cvZero( im );

		printf("Fitting PROSAC\n");
		CvCircle32f c = cvFitCircleRobust1( &data, &weights, &inliers, 200000, .01);
		printf("Fitting RANSAC\n");
		CvCircle32f c2 = cvFitCircleRobust( &data, 200000, .01);
		printf("%f %f %f\n", c.c.x, c.c.y, c.r);

		for(int i=0; i<data.rows; i++){
			if(inliers[i][0]>0){
				cvCircle(im, cvPoint( data[i][0], data[i][1] ), 3, CV_RGB(255,0,0));
			}
			else{
				cvCircle(im, cvPoint( data[i][0], data[i][1] ), 3, CV_RGB(255,255,255));
			}

		}

		cvCircle(im, cvPointFrom32f( c.c ), (int) c.r, CV_RGB(0, 255, 0) );
		cvCircle(im, cvPointFrom32f( c2.c ), (int) c2.r, CV_RGB(0, 255, 255) );
		cvShowImage("win", im);
		cvWaitKey(-1);
	}
}
