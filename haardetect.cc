#include <cv/haar_detector.hpp>

int main(int argc, char ** argv){
	HaarDetector haar(argv[1]);
	cv::Image<uchar, 1> im;
	cv::Image<uchar, 3> im3;
	cvNamedWindow("win", 0);
	for(int i=2; i<argc; i++){
		assert(im.open(argv[i]));
		im3 = im;

		haar.detect(im);
	
		haar.draw(&im3, CV_RGB(0,0,255));

		im3.show("win");
		cvWaitKey(-1);
	}

}
