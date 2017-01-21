#include <unistd.h>
#include "im_mmap.h"
#include "cv.h"
#include "highgui.h"

int main(int argc, char ** argv){
	if(argc<3) return -1;
	int fd = im_mmap_open(argv[1], O_CREAT);
	if(fd<0) return -1;
	IplImage * im = cvLoadImage(argv[2]);
	write(fd, im, sizeof(IplImage));
	write(fd, im->imageData, im->imageSize);
	close(fd);
}
