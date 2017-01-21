#include "stereo_camera.hpp"

#define FLEA_GREY 0x00b09d01004ff57f
#define FLEA_HICOLOR 0x00b09d0100515c58

int main(int argc, char ** argv){
	StereoCamera cam(FLEA_HICOLOR, FLEA_GREY);
	cvNamedWindow("win0", 0);
	cvNamedWindow("win1", 0);
	cv::Image<uchar> flipped;
	while(1){
		cam.lock();
		cam.color().resize(640,480).show("win0");
		flipped.realloc( cam.right() );
		cvFlip( &cam.right(), &flipped, 1 );
		flipped.resize(640,480).show("win1");
		if((char)cvWaitKey(10)=='q') break;
	}
}
#if 0
#include <libdc1394++/DC1394Camera.hpp>
#include <cv/image.h>

int main(int argc, char ** argv){
	DC1394::Camera cam0(0,0);
	DC1394::Camera cam1(0,1);
	IplImage header;

	assert(cam0.setup(-1, FORMAT_SVGA_NONCOMPRESSED_1, MODE_1024x768_MONO, SPEED_400, FRAMERATE_15));
	assert(cam1.setup(-1, FORMAT_SVGA_NONCOMPRESSED_1, MODE_1024x768_MONO, SPEED_400, FRAMERATE_15));
	assert(cam0.isValid());
	assert(cam0.start());
	assert(cam1.isValid());
	assert(cam1.start());

	cvNamedWindow("win0", 0);
	cvNamedWindow("win1", 0);
	
	cvInitImageHeader(&header, cvSize(cam0.getWidth(), cam0.getHeight()), 8, 1);

	while(1){
		assert( cam0.capture() );
		assert( cam1.capture() );
		header.imageData = header.imageDataOrigin = (char *) cam0.getBuffer();
		cvShowImage("win0", &header);
		header.imageData = header.imageDataOrigin = (char *) cam1.getBuffer();
		cvShowImage("win1", &header);
		if((char)cvWaitKey(10)=='q') break;
	}
}
#endif
