#include <cv/image.h>
#include <cxcore.h>
#include <cv.h>
#include <highgui.h>

int main(int argc, char ** argv){
	cv::Image<float> im;
	cv::Image<float> rot;
	im.open( argv[1] );
	rot.realloc( im );

	float theta_step = M_PI/32;
	int dft_M = cvGetOptimalDFTSize( im.height - 1 );
	int dft_N = cvGetOptimalDFTSize( im.width - 1 );

	cv::Image<float> dft_A( dft_N, dft_M );
	CvMat tmp;
	cvNamedWindow("win", 0);

	int thresh=0;
	cvCreateTrackbar("Threshold", "win", &thresh, 100, NULL); 
	while(1){

		// copy A to dft_A and pad dft_A with zeros
		cvGetSubRect( &dft_A, &tmp, cvRect(0,0,im.width, im.height));
		cvCopy(&im, &tmp);

		cvGetSubRect( &dft_A, &tmp, cvRect(im.width, 0 , dft_A.width-im.width, dft_A.height) );
		cvZero( &tmp );
		// no need to pad bottom part of dft_A with zeros because of
		// use nonzero_rows parameter in cvDFT() call below
		cvDFT( &dft_A, &dft_A, CV_DXT_FORWARD + CV_DXT_SCALE, im.height );

		double min, max;
		cvMinMaxLoc( &dft_A, &min, &max );
		printf("minmax = %f %f\n", min, max);

		cvThreshold( &dft_A, &dft_A, thresh*(max-min)*.01+min, 0, CV_THRESH_TOZERO );

		cvDFT( &dft_A, &dft_A, CV_DXT_INV_SCALE, im.height );
		
		//cvLog(&dft_A, &dft_A);
		cvGetSubRect( &dft_A, &tmp, cvRect(0,0,im.width, im.height));
		cvCopy(&tmp, &rot);

		rot.imagesc("win");
		cvWaitKey(10);
	}
}
