#include <cv/image.h>
#include <cv/gaussian_scale_space.hpp>

template< typename data_t >
class GeometricBlur : public cv::GaussianPyramid<data_t> {
public:
	GeometricBlur() : cv::GaussianPyramid<data_t>(16, 5) {}

	void blur(int x, int y, cv::Image<data_t,1> & res, float alpha){
		cv::Image<data_t, 1> tmp = this->getImage(0,0); 
		cv::Image<uchar, 1> mask = tmp;
		res = tmp;
		
		// initialize
		cvResize(&(this->getImage(this->getNumScales()-1, this->getNumSubScales()-1)), &tmp);
		cvCopy(&tmp, &res);

		for(int i=this->getNumScales()-1; i>=0; i--){
			for(int j=this->getNumSubScales()-1; j>=0; j--){
				float sigma = this->getSigma(i,j); 
				// mask is a circle with radius sigma * alpha	
				cvZero(&mask);
				cvResize(&(this->getImage(i,j)), &tmp);
				printf("Radius = %f\n", sigma*alpha);
				cvCircle(&mask, cvPoint(x,y), sigma*alpha, cvScalarAll(255), -1);
				cvCopy(&tmp, &res, &mask);
			}
		}
	}
	void blur2(int x, int y, cv::Image<data_t,1> & res, float alpha){
		cv::Image<data_t, 1> tmp = (*this)[0]; 
		cv::Image<uchar, 1> mask = tmp;
		res = tmp;
		cvSet(&mask, cvScalarAll(255));
		for(float r=1; r<15; r++){
			float sigma = alpha*r;
			cvSmooth(&tmp, &tmp, CV_GAUSSIAN, 0, 0 ,sigma);
			cvCircle(&mask, cvPoint(x,y), r, cvScalarAll(0), -1);
			cvCopy(&tmp, &res, &mask);
		}
	}
};
int mouseX, mouseY;
void onMouse( int event, int x, int y, int flags , void * arg){
	mouseX = x;
	mouseY = y;
}


int main(int argc, char ** argv){
	cv::Image<uchar,1> im;
	GeometricBlur <uchar> pyr;
	mouseX = 0;
	mouseY = 0;

	assert(im.open(argv[1]));

	pyr.calculate(im);

	cvNamedWindow("win", 0);
	cvSetMouseCallback("win", onMouse, &pyr);

	while(1){
		pyr.blur( mouseX, mouseY, im, 2);
		im.show("win");
		if((uchar)cvWaitKey(100)=='q') break;
	}
}
