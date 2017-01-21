#include <cv/circle_detector.hpp>

#define L2C(c1, c2) sqrt( (c2.c.x-c1.c.x)*(c2.c.x-c1.c.x) + (c2.c.y-c1.c.y)*(c2.c.y-c1.c.y) +(c2.r-c1.r)*(c2.r-c1.r)) 

void ground_truth( cv::Image<float> & im, CvCircle32f * circles, int num, int zoom=5){
	cv::Image<float> zim(im.width*zoom, im.height*zoom);
	cvSet( &zim, cvScalarAll(255) );
	for(int i=0; i<num; i++){
		CvPoint p = cvPoint( drand48()*zim.width, drand48()*zim.height );
		int r = (int) 10 + (drand48()*( MIN(zim.width,zim.height)/10 ));
		cvCircle( &zim, p, r, cvScalarAll(0), -1, CV_AA);

		circles[i].c.x = p.x*(1.0/zoom);
		circles[i].c.y = p.y*(1.0/zoom);
		circles[i].r = r*(1.0/zoom);
	}
	for(int j=0; j<zim.height; j++){
		for(int i=0; i<zim.width; i++){
			float v = zim[j][i] + (drand48()*510 - 255);
			if(v<0) v = 0;
			if(v>255) v = 255;
			zim[j][i] = v;
		}
	}
	
	cvSmooth( &zim, &zim, CV_GAUSSIAN, 11, 11);
	cvResize( &zim, &im, CV_INTER_CUBIC );
}
int min_l2( CvCircle32f * circles, int ncircles, CvCircle32f c){
	int minidx=0;
	float mindist=L2C(circles[0], c); 
	for(int i=1; i<ncircles; i++){
		float d = L2C(circles[i], c);
		if(d<mindist){
			mindist=d;
			minidx=i;
		}
	}
	return minidx;
}
	
int main(int argc, char ** argv){
	CircleDetector circledetector;
	cv::Image<float> im(640,480);
	cv::Image<float> Iy(640,480);
	cv::Image<float> Ix(640,480);
	cv::Image<uchar,3> im3;
	CvCircle32f circles[5];
	int ncircles = 5;

	int n=0;
	double sum1=0, sum2=0, sum3=0, sum4=0;
	double ssq1=0, ssq2=0, ssq3=0, ssq4=0;
	cvNamedWindow("win", 0);
	while(1){
		ground_truth( im, circles, ncircles, 5 );
		cvSobel(&im, &Ix, 1, 0, CV_SCHARR );
		cvSobel(&im, &Iy, 0, 1, CV_SCHARR );

		circledetector.detect( &im, 128 );
		im3 = im;

		//circledetector.getThreshImage().show("win");
		//cvWaitKey(100);
		for(int i=0; i<circledetector.size(); i++){
			CvCircle32f c = circledetector[i];
			CvCircle32f c1,c2;
			int idx = min_l2( circles, ncircles, c);
			float d4, d1, d2, d = L2C(circles[idx], c);
			if(d > 5) continue;
			//printf("%f, %f, %f\n", circles[idx].c.x, circles[idx].c.y, circles[idx].r);
			//printf("%f, %f, %f : %f\n", c.c.x, c.c.y, c.r, d );
			cvCircle(&im3, cvPointFrom32f( c.c ), (int) c.r, CV_RGB(0,255,0), 1, CV_AA );
			
			// gradient optimization
			c1 = CircleDetector::refine( c, Ix, Iy );
			cvCircle(&im3, cvPointFrom32f( c1.c ), (int) c1.r, CV_RGB(0,0,255), 1, CV_AA );
			d4 =  L2C(circles[idx], c1); // + pow(circles[idx].r-c.r, 2);

			c2 = CircleDetector::refine( c, im );
			cvCircle(&im3, cvPointFrom32f( c2.c ), (int) c2.r, CV_RGB(255,0,0), 1, CV_AA );
			d1 = L2C(circles[idx], c2); // + pow(circles[idx].r-c.r, 2);
			//printf("%f, %f, %f : %f\n", c.c.x, c.c.y, c.r, d1 );

			
			c = CircleDetector::subpixel_refine( c1, Ix, Iy);
			cvCircle(&im3, cvPointFrom32f( c.c ), (int) c.r, CV_RGB(255,0,0), 1, CV_AA );
			d2 = L2C(circles[idx], c); // + pow(circles[idx].r-c.r, 2);
			if(d4!=d2){
				printf("%f, %f, %f : %f %f %f\n", c.c.x, c.c.y, c.r, d, d1, d2 );
				printf("baseline avg: %f (%f)\n", sum3/n, sqrt( ssq3/n - sum3*sum3/(n*n) ));
				printf("refinement avg: %f (%f)\n", sum1/n, sqrt( ssq1/n - sum1*sum1/(n*n) ));
				printf("gradient refinement avg: %f (%f)\n", sum4/n, sqrt( ssq4/n - sum4*sum4/(n*n) ));
				printf("subpixel avg: %f (%f)\n", sum2/n, sqrt( ssq2/n - sum2*sum2/(n*n) ));
				printf("subpixel improvement: %f\n", (sum4-sum2)/n );
			}
			sum1+=d1;
			sum2+=d2;
			sum3+=d;
			sum4+=d4;
			ssq1+=(d1*d1);
			ssq2+=(d2*d2);
			ssq3+=(d*d);
			ssq4+=(d4*d4);
			n++;

		}
		im3.show("win");
		cvWaitKey(100);
	}
}
