#include <stdio.h>

/* According to POSIX 1003.1-2001 */
#include <sys/select.h>

/* According to earlier standards */
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
										  
#include "cv.h"
#include "highgui.h"

int main(int argc, char ** argv){
	IplImage im;
	// select on stdin
	fd_set rfds;
	struct timeval tv;
	int retval;
	im.imageData = (char *) malloc(16);
	cvNamedWindow("win", 0);
	printf("piped_viewer started %s\n", argv[1]);
	do {
		/* Watch stdin (fd 0) to see when it has input. */
		FD_ZERO(&rfds);
		FD_SET(0, &rfds);
		
		tv.tv_sec = 0;
		tv.tv_usec = 1;

		retval = select(1, &rfds, NULL, NULL, &tv);

		if (retval == -1)
			perror("select()");
		else if(retval!=0){
			if(fread(&im, sizeof(IplImage), 1, stdin)!=1){
				//perror("fread()");
				return -1;
			}
			cvInitImageHeader(&im, cvSize(im.width, im.height), im.depth, im.nChannels);
			cvCreateData(&im);
			if(im.imageSize<=0){
				fprintf(stderr, "WARNING: zero length image passed");
				continue;
			}
			if(fread(im.imageData, 1, im.imageSize, stdin) != im.imageSize){
				//perror("fread()");
				return -1;
			}
			cvShowImage("win", &im);
		}	
	}while(cvWaitKey(10)!='q');
}
