#include <cv/image.h>
#include <cv/contour_detector.hpp>
#include <highgui.h>
#include <iostream>

int main(int argc, char ** argv){
	cv::Image<uchar> im,count,min,max;
	cv::Image<uchar, 3 > im3, cp3;
	im3.open(argv[1]);
	cv::Image<float> f1,f2,d1,d2,d3;
	int gwidth = 3;
	double sigma = (3*.5-1.0)*.3 + .8;

	//im3 = im3.scale(100,100);
	im = im3;
	f1 = im;
	count = im;
	min = im;
	max = im;
	d3 = f1;
	d2 = d3;
	d1 = d2;
	
	cvNamedWindow("w", 0);
	cvNamedWindow("d1", 0);
	cvNamedWindow("d2", 0);
	cvNamedWindow("d3", 0);
	ContourDetector contours;
	cp3=im3;
	for(int i=0; i<100000; i++){
		d3 = d2;
		d2 = d1;
		f2 = f1;
		
		cvSmooth(&f1, &f1, CV_GAUSSIAN, gwidth, gwidth);
		d1 = f1 - f2;
		cvAbs(&d1, &d3);

		d3.normalize();
		contours.detect(&d3, 64.0, false);
		//cp3=im3;
		contours.draw(&cp3, CV_RGB(i%255, i/255, 255));

		cp3.show("w");
		d1.imagesc("d1");
		d2.imagesc("d2");
		d3.imagesc("d3");
		
		std::cout<<i<<std::endl;
		if(cvWaitKey(10)=='q') break;
	}	
	cvWaitKey(-1);
}
