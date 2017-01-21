#include <stdio.h>
#include <cxcore.h>
#include <cv.h>
#include <highgui.h>

#define DEFAULT_CODEC "DIVX"
#define DEFAULT_FPS 30
#define DEFAULT_BPS 1000

void usage(char ** argv ){
	fprintf(stderr, "%s <output file> [ codec ] [ framerate ] [ bitrate ]\n", argv[0]);
}
int main(int argc, char ** argv){
	const char * codec = DEFAULT_CODEC;
	double fps = DEFAULT_FPS;
	int bps = DEFAULT_BPS;
	CvSize sz = cvSize(640, 480);

	if(argc<2){
		usage(argv);
		return -1;
	}
	if(argc>4) bps = atoi(argv[4]);
	if(argc>3) fps = atof(argv[3]);
	if(argc>2) codec = argv[2];

	CvCapture * cap = cvCaptureFromCAM(0);
	assert(cap);
	IplImage * im = cvQueryFrame( cap );
	assert(im);
	IplImage * resize_im = cvCreateImage( sz, 8, 3 );
	CvVideoWriter * videowriter = cvCreateVideoWriter( argv[1], CV_FOURCC(codec[0], codec[1], codec[2], codec[3]), fps, sz);
	cvNamedWindow("win", 0);
	while(1){
		im = cvQueryFrame( cap );
		cvResize(im, resize_im);
		cvShowImage("win", resize_im);
		cvWriteFrame(videowriter, resize_im);
		if((char)cvWaitKey(10)=='q') break;
	}
	cvReleaseVideoWriter(&videowriter);
}
