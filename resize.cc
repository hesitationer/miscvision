#include <cv/image.h>
void resize(cv::Image<float,1> & src, cv::Image<float, 1> & dst, int width, int height){
	float vstep=(double)src.width/(width);
	float hstep=(double)src.height/(height);
	dst.reshape(width,height);
	for(int j=0; j<height; j++){
		float * ptr = dst[j];
		for(int i=0; i<width; i++){
			ptr[i] = src.getSubPix(i*vstep, j*hstep);
		}
	}
}
void decimate(cv::Image<float,1> &src, cv::Image<float,1> & dst, int factor){
	int width=dst.width;
	int height=dst.height;
	for(int j=0; j<height; j++){
		for(int i=0; i<width; i++){
			dst[j][i]=src[j*factor][i*factor];
		}
	}
}

int main(int argc, char ** argv){
	cv::Image<float, 1> orig, r1, r2;
	assert(orig.open(argv[1]));
	r1.realloc(orig.width/2, orig.height/2);
	r2.realloc(orig.width/2, orig.height/2);
	
	cvResize(&orig, &r1, CV_INTER_NN);
	decimate(orig, r2, 2);

	cvNamedWindow("win", 0);
	r1.show("win"); cvWaitKey(-1);
	r2.show("win"); cvWaitKey(-1);

	r1 = r1-r2;
	printf("%d\n", cvCountNonZero( &r1 ) );
	r1.show("win"); cvWaitKey(-1);
}


