#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "cv.h"
#include "highgui.h"

int main(int argc, char ** argv) {
	FILE *command,*data;
	char *temp_name = NULL;
	double a,b;
	int i;
	char command_string[256];
	IplImage * im = cvLoadImage(argv[2]);
	if(im==NULL) return -1;
	sprintf(command_string, "/apt1/cs/student/rstancha/projects/research/cv/piped_viewer %s", argv[1]);
	command = popen(command_string,"w");
	printf("writing image header\n");
	fflush(command);
	fwrite(im, sizeof(IplImage),  1, command);
	fflush(command);
	printf("writing %d x %d image = %d\n", im->width, im->height, im->imageSize);
	fwrite(im->imageData, 1,  im->imageSize, command);
	fflush(command);
	//fprintf(command, "\n");
	fprintf(stderr, "press enter to continue..."); fflush(stderr);
	pclose(command);
	return 0;
}
