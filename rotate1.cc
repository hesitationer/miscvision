#include <stdio.h>
#include <highgui.h>
#include <cv/cvext.h>

CvPoint points[2];
int npoints=0;

void cb( int event, int x, int y, int flags, void * param ){
	if(event==CV_EVENT_LBUTTONUP){
		if(npoints<2){
			points[npoints] = cvPoint(x,y);
			npoints++;
		}
	}
}
#define DIST(p,q) (sqrt( (p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y)))

int main(int argc, char ** argv){
	cvNamedWindow("image", 0);
	cvNamedWindow("rotated", 0);
	IplImage * im = cvLoadImage( argv[1], 1 );
	IplImage * clone  = cvCloneImage( im );
	IplImage * rotated = cvCreateImage(cvSize(200,200), 8, 3); 

	cvSetMouseCallback("image", cb);
	while(1){
		if(npoints==1){
			cvCopy(im, clone);
			cvCircle( clone, points[0], 1, CV_RGB(0,255,0) );
		}
		if(npoints==2){
			float d = DIST( points[0], points[1] );
			float angle = atan2( points[1].y-points[0].y, points[1].x-points[0].x );
			float scale = rotated->width*.5/d;
			printf("scale=%f angle=%f\n", scale, angle/M_PI*180);
			cvRotateZoomImage( im, rotated, cvPointTo32f(points[0]), angle, scale );
			cvShowImage("rotated", rotated);
			cvCircle( clone, points[1], 1, CV_RGB(0,255,0) );
			cvLine( clone, points[0], points[1], CV_RGB(0,255,0) );
			npoints=0;
		}
		cvShowImage("image", clone);
		if((char)cvWaitKey(100)=='q') break;
	}
}
