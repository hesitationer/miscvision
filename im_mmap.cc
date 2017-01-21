#include <unistd.h>
#include "im_mmap.h"
#include "highgui.h"

int main(int argc, char ** argv){
	if(argc<2){
		return -1;
	}
	cvNamedWindow(argv[1], 0);
	int fd = -1;
	while(1){
		fd = im_mmap_open(argv[1], 0);
		if(fd>=0) break;
		usleep(50000);
	}
	do{
		IplImage * im = im_mmap_read(fd);
		cvShowImage(argv[1], im);
		cvWaitKey(10);
	}while(cvGetWindowHandle(argv[1])!=NULL);
}
