#include <cv/image.h>
#include <cv/matrix.hh>
#include <cv/cvext.h>
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
	cvNamedWindow("fft", 0);
	cvNamedWindow("win", 0);
	for(float theta=0; theta<M_PI*2; theta+=theta_step){

		cvRotateZoomImage(&im, &rot, cvPoint2D32f( im.width/2, im.height/2 ), theta, 1);
		
		// copy A to dft_A and pad dft_A with zeros
		cvGetSubRect( &dft_A, &tmp, cvRect(0,0,im.width, im.height));
		cvCopy(&rot, &tmp);
		cvGetSubRect( &dft_A, &tmp, cvRect(im.width, 0 , dft_A.width-im.width, dft_A.height) );
		cvZero( &tmp );
		// no need to pad bottom part of dft_A with zeros because of
		// use nonzero_rows parameter in cvDFT() call below

		cvDFT( &dft_A, &dft_A, CV_DXT_FORWARD, im.height );


		//cvLog(&dft_A, &dft_A);

		rot.imagesc("win");
		dft_A.imagesc("fft");
		cvWaitKey(-1);
	}
}
