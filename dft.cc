#include <cxcore.h>
#include <cv.h>
#include <cv/cvext.h>
#include <highgui.h>


int main(int argc, char ** argv){
	cvRedirectError( cvSegReport );
	IplImage * im = cvLoadImage( argv[1], CV_LOAD_IMAGE_GRAYSCALE);

	IplImage * realInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 1);
	IplImage * imaginaryInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 1);
	IplImage * complexInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 2);

	cvScale(im, realInput);
	cvZero(imaginaryInput);
	cvMerge(realInput, imaginaryInput, NULL, NULL, complexInput);

	int dft_M = cvGetOptimalDFTSize( im->height - 1 );
	int dft_N = cvGetOptimalDFTSize( im->width - 1 );

	CvMat* dft_A = cvCreateMat( dft_M, dft_N, CV_64FC2 );
	IplImage * image_Re = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);
	IplImage * image_Im = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);
	CvMat tmp;

	// copy A to dft_A and pad dft_A with zeros
	cvGetSubRect( dft_A, &tmp, cvRect(0,0, im->height, im->width));
	cvCopy( complexInput, &tmp );
	cvGetSubRect( dft_A, &tmp, cvRect(im->width,0, dft_A->cols - im->width, im->height));
	cvZero( &tmp );

	// no need to pad bottom part of dft_A with zeros because of
	// use nonzero_rows parameter in cvDFT() call below

	cvDFT( dft_A, dft_A, CV_DXT_FORWARD, complexInput->height );

	cvNamedWindow("win", 0);
	cvNamedWindow("magnitude", 0);
	cvNamedWindow("phase", 0);
	cvShowImage("win", im);

	// Split Fourier in real and imaginary parts
	cvSplit( dft_A, image_Re, image_Im, 0, 0 );

	// Compute the magnitude of the spectrum Mag = sqrt(Re^2 + Im^2)
	cvPow( image_Re, image_Re, 2.0);
	cvPow( image_Im, image_Im, 2.0);
	cvAdd( image_Re, image_Im, image_Re);
	cvPow( image_Re, image_Re, 0.5 );

	// Compute log(1 + Mag)
	cvAddS( image_Re, cvScalar(1.0), image_Re ); // 1 + Mag
	cvLog( image_Re, image_Re ); // log(1 + Mag)


	// Rearrange the quadrants of Fourier image so that the origin is at
	// the image center
	cvShiftDFT( image_Re, image_Re );

	// Localize( image_Re, &minVal, &maxVal );	
	double min,max;
	cvMinMaxLoc(image_Re, &min, &max);
	cvScale(image_Re, image_Re, 1.0/(max-min), 1.0*(-min)/(max-min));
	cvShowImage("magnitude", image_Re);

	cvWaitKey(-1);
}
