#include <cv/cvext.h>
#include <cv/image.h>

void draw(IplImage * im, cv::Image<uchar,1> &mask, CvScalar color, int radius){
	for(int j=0; j<mask.height; j++){
		for(int i=0; i<mask.width; i++){
			if(mask[j][i]==0) continue;
			cvCircle(im, cvPoint(i*im->width/mask.width,j*im->height/mask.height), radius, color, 1, CV_AA);
		}
	}
}

int main(int argc, char ** argv){
	cv::Image<float, 1> im, im2;
	cv::Image<uchar, 3> im3;
	cv::Image<uchar, 1> min,max;
	cv::Image<float, 1> *g0, *g1, *g2, *g3;
	cv::Image<float, 1> *d0, *d1, *d2, *tmp;

	im.open(argv[1]);
	min.reshape(im);
	max.reshape(im);
	g0 = new cv::Image<float, 1>( im );
	g1 = new cv::Image<float, 1>( im );
	g2 = new cv::Image<float, 1>( im );
	g3 = new cv::Image<float, 1>( im );
	d0 = new cv::Image<float, 1>( im );
	d1 = new cv::Image<float, 1>( im );
	d2 = new cv::Image<float, 1>( im );
	float sigma = ((5/2-1)*0.3 + 0.8)/2;  // pyr down is 5x5 gaussian...
	// want pyr down to be 2*sigma... hence the /2 term
	float sigmastep = exp2(1./10);
	float esigma = 1;

	// construct initial scales
	cvSmooth(&im, g0, CV_GAUSSIAN, 0, 0, sigma); esigma *= sigma*sigmastep;
	cvSmooth(&im, g1, CV_GAUSSIAN, 0, 0, esigma); esigma *= sigmastep;
	cvSub(g0, g1, d0);
	cvSmooth(&im, g2, CV_GAUSSIAN, 0, 0, esigma); esigma *= sigmastep;
	cvSub(g1, g2, d1);
	cvNamedWindow("win", 1);
	while(im.width >= 16 && im.height >= 16){
		im3.open(argv[1]);
		cvSmooth(&im, g3, CV_GAUSSIAN, 0, 0, esigma); esigma *= sigmastep;
		cvSub(g2, g3, d2);
		cvIsLocalMin2(*d1, *d1, min, false, 3);
		cvIsLocalMin2(*d1, *d0, min, true, 3, &min);
		cvIsLocalMin2(*d1, *d2, min, true, 3, &min);
		cvIsLocalMax2(*d1, *d1, max, false, 3);
		cvIsLocalMax2(*d1, *d0, max, true, 3, &max);
		cvIsLocalMax2(*d1, *d2, max, true, 3, &max);
		
//		cvOr(&min, &max, &max);

		int radius = ((((3.0 * esigma/(sigmastep*sigmastep)) * 3.0) / 2.0) + 0.5)/2;
		draw(&im3, min, CV_RGB(0,255,0), radius);
		draw(&im3, max, CV_RGB(0,0, 255), radius);
		im3.show("win");
		cvWaitKey(-1);
		printf("effective sigma = %f\n", esigma);
		if(esigma/(sigmastep*sigmastep) >= 2*sigma){
			im2 = im;
			im.reshape(im.width/2, im.height/2);
			cvPyrDown(&im2, &im);
			g0->reshape(im);
			g1->reshape(im);
			g2->reshape(im);
			g3->reshape(im);
			d0->reshape(im);
			d1->reshape(im);
			d2->reshape(im);
			min.reshape(im);
			max.reshape(im);
			sigma*=2;
			(*g0) = im;	esigma = sigma*sigmastep;
			cvSmooth(&im, g1, CV_GAUSSIAN, 0, 0, esigma); esigma *= sigmastep;
			cvSub(g0, g1, d0);
			cvSmooth(&im, g2, CV_GAUSSIAN, 0, 0, esigma); esigma *= sigmastep;
			cvSub(g1, g2, d1);
			cvWaitKey(-1);
		}


		tmp=g0;
		g0=g1;
		g1=g2;
		g2=g3;
		g3=tmp;
		tmp=d0;
		d0=d1;
		d1=d2;
		d2=tmp;
	}
}
