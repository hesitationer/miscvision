#include <cv/image.h>
#include <cv/gaussian_scale_space.hpp>

int main(int argc, char ** argv){
	cv::Image<uchar,1> im;
	cv::Image<float,1> Ix, Iy;
	cv::GaussianPyramid <float> pyr(16,3);

	assert(im.open(argv[1]));

	pyr.calculate(im);

	cvNamedWindow("win", 0);
	for(int i=0; i<pyr.getNumScales(); i++){
		for(int j=0; j<pyr.getNumSubScales(); j++){
			Ix.reshape(pyr.getImage(i,j));
			Iy.reshape(pyr.getImage(i,j));
			cvSobel( &(pyr.getImage(i,j)), &Ix, 1, 0, CV_SCHARR);
			cvSobel( &(pyr.getImage(i,j)), &Iy, 1, 0, CV_SCHARR);
			Ix = (Ix*Ix) + (Iy*Iy);
			Ix.imagesc( "win" );
			if((uchar)cvWaitKey(-1)=='q') break;
		}
	}
}
