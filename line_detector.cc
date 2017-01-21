#include <cv/image.h>
#include <cv/line_detector.hpp>

int main(int argc, char ** argv){

	cv::Image<uchar, 1> im, dst;
	cv::Image<uchar, 3> im3;
	LineDetector linedetector;
	for(int i=1; i<argc; i++){
		assert( im.open( argv[i] ) );
		im3 = im;
		dst.realloc( im );
	
		cvNamedWindow("win", 0);
		cvCanny( &im, &dst, 50, 200, 3 );
		dst.show("win"); cvWaitKey();
		linedetector.detect( &dst);
		linedetector.draw( &im3, CV_RGB(0,255,0) );
		im3.show("win");
		cvWaitKey();
	}
}

