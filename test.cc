#include <cv.h>  // basic cv includes
#include <highgui.h> // stuff for gui's, opening images

int main(int argc, char ** argv){

	// load image given as argument
	IplImage * im = cvLoadImage(argv[1]);

	// Create a window
	cvNamedWindow("win", 0);
	
	// Show the image
	cvShowImage("win", im);

	// Wait for response indefinitely
	cvWaitKey(-1);
}
