#include<cv/image.h>
#include <cv/cvext.h>
int main(int argc, char ** argv){
	cv::Image<uchar, 1> im;
	cv::Image<float, 1> harris;
	assert(im.open(argv[1]));
	harris = im;
	cvCornerHarris( &im, &harris, 11, 3, 0.04);

	cvIsLocalMax( &harris, &im, 5 );
	cvNamedWindow("win", 1);
	im.show("win");//.imagesc("win", -0.0001, 0.0001);
	cvWaitKey(-1);

	im.open(argv[1]);
	cvPreCornerDetect( &im, &harris);
	cvIsLocalMax( &harris, &im, 5);
	im.show("win");//.imagesc("win", -0.0001, 0.0001);
	cvWaitKey(-1);

	im.open(argv[1]);
	cvCornerMinEigenVal( &im, &harris, 5 );
	cvIsLocalMax( &harris, &im );
	im.show("win");
	cvWaitKey(-1);


}
