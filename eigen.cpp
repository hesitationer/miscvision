#include <stdio.h>
#include <cv/timer.h>
#include <cv/cvext.h>

int main(){
	CvMat * EV1 = cvCreateMat(6, 256*256, CV_32F);
	CvMat * EV2 = cvCreateMat(500, 256*256, CV_32F);
	CvMat * coef = cvCreateMat(3, 1, CV_32F);
	CvMat * data = cvCreateMat(256*256, 1, CV_32F);
	CvMat * mean = cvCreateMat(256*256, 1, CV_32F);
	cvSetIdentity( EV1 );
	cvSetIdentity( EV2 );

	Timer timer;
	timer_start(&timer);
	for(int i=0; i<1000; i++){
		cvEigenBackProject(coef, mean, EV1, data);
	}
	printf("%fs\n", timer_elapsed(&timer)/1000 );
	timer_start(&timer);
	for(int i=0; i<1000; i++){
		cvEigenBackProject(coef, mean, EV2, data);
	}
	printf("%fs\n", timer_elapsed(&timer)/1000 );
}
