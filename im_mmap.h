#ifndef __IM_MMAP_H
#define __IM_MMAP_H
#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "cv.h"

int im_mmap_open(char * name, int flags){
	char fname[256];
	sprintf(fname, "/tmp/im_mmap_%s", name);
	return open(fname, O_RDWR | flags,  S_IRWXU);
}

int im_mmap_write(int fd, IplImage * im){
	lseek(fd, sizeof(IplImage)+im->imageSize, SEEK_SET);
	char * ptr = (char *) mmap(0, sizeof(IplImage)+im->imageSize, PROT_WRITE, MAP_SHARED, fd, 0);
	if(!ptr) return -1;
	memcpy(ptr, im, sizeof(IplImage));
	memcpy(ptr+sizeof(IplImage), im->imageData, im->imageSize);
	return 0;
}

IplImage * im_mmap_read(int fd){
	IplImage * im = (IplImage *) mmap(0, sizeof(IplImage), PROT_READ, MAP_SHARED, fd, 0);
	if(!im) return NULL;
	im->imageData = (char *) mmap(0, im->imageSize, PROT_READ, MAP_SHARED, fd, sizeof(IplImage));
	return im;
}

#endif
