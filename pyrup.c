#include <stdio.h>
#include <cv.h>
#include <highgui.h>

int main(int argc, char ** argv){
	IplImage * im1 = cvLoadImage(argv[1],0);
	IplImage * im2 = cvCreateImage(cvSize(im1->width*2, im1->height*2), im1->depth, im1->nChannels);
	cvPyrUp(im1,im2);
	cvSaveImage("pyrup.bmp", im2);
	return 0;
}
