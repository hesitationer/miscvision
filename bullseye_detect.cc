#include <cv/image.h>
#include <cv/matlab_compat.h>
#include <cv/bullseye_detector2.hh>
int main(int argc, char ** argv){
	BullsEyeDetector bd;
	matlab_compat::figure(1);
	cv::Image<uchar> im;
	cv::Image<uchar,3> im3;
	assert( im.open(argv[1]) );

	
	cv::Image<uchar> im1;
	im1.reshape(im);
	

	bd.detect(&im);

	im3 = im;
	bd.draw(&im3, CV_RGB(0,255,0));

	matlab_compat::imshow(&im);
	matlab_compat::drawnow();
	matlab_compat::ginput();	
	matlab_compat::imshow(&im3);
	matlab_compat::drawnow();
	matlab_compat::ginput();	
}
