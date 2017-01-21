#include <cv/corner_detector.hpp>

int main(int argc, char ** argv){
	int method=0;
	int thresh=0;
	cv::Image<uchar, 1> im;
	im.open(argv[1]);
	cv::Image<uchar, 3> im3;
	CornerDetector detector;
	std::vector< CvPoint > corners;
	cvNamedWindow("win", 0);
	int win = 3;
	cvCreateTrackbar("Window", "win", &win, 100, NULL);
	cvCreateTrackbar("Threshold", "win", &thresh, 255, NULL);
	cvCreateTrackbar("Method", "win", &method, 2, NULL);

	while(1){
		if(win<3) win=3;
		if(win%2==0) win++;
		detector.detect(&im, corners, thresh+1, win, method);
		im3 = im;
		for(int i=0; i<corners.size(); i++){
			cvCircle(&im3, corners[i], 1, CV_RGB(0, 255, 0), 1, CV_AA );
		}
		im3.show("win");
		if((char)cvWaitKey()=='q') break;
	}

}
