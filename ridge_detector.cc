#include <cv/image.h>
#include <cv/line_detector.hpp>
#include <cv/image_pyramid.hh>
#include <cv/ridge_detector.hpp>

int gx,gy;
void cb(int event, int x, int y, int flags, void * userdata){
	gx = x;
	gy = y;
}
void draw( IplImage * im, RidgeDetector & r, float scale, CvScalar color){
	float rad = 2;
	for(int i=0; i<r.size(); i++){
		CvPoint p;
		CvPoint p1, p2;
		float theta;
		r.get(i, &p, &theta);
		p1 = cvPoint( scale*(p.x - cos(theta)*rad), scale*(p.y - sin(theta)*rad) );
		p2 = cvPoint( scale*(p.x + cos(theta)*rad), scale*(p.y + sin(theta)*rad) );
		p1.x -= scale*sin(theta)/2;
		p1.y += scale*cos(theta)/2;
		p2.x -= scale*sin(theta)/2;
		p2.y += scale*cos(theta)/2;
		cvLine(im, p1, p2, color, 1, CV_AA); 
		p1.x += scale*sin(theta);
		p1.y -= scale*cos(theta);
		p2.x += scale*sin(theta);
		p2.y -= scale*cos(theta);
		cvLine(im, p1, p2, color, 1, CV_AA); 
	}
}
void draw_clusters( IplImage * im, std::vector<CvPoint2D32f> & lines, float scale, CvScalar color ){
	CvSize sz = cvGetSize(im);
	for(int i=0; i<lines.size(); i++){
		float rho = lines[i].x*scale;
		float theta = lines[i].y;
		CvPoint2D32f p1,p2;
		if(LineDetector::fromHoughLine( rho-scale/2, theta, sz, &p1, &p2 )){
			cvLine(im, cvPointFrom32f( p1 ), cvPointFrom32f( p2 ), color, 1, CV_AA);
		}
		if(	LineDetector::fromHoughLine( rho+scale/2, theta, sz, &p1, &p2 )){
			cvLine(im, cvPointFrom32f( p1 ), cvPointFrom32f( p2 ), color, 1, CV_AA);
		}
	}
}
IplImage * gimage;
float gscale=1.0;
void cluster( RidgeDetector & r, int rhosize, std::vector<CvPoint2D32f> & lines){
	int thetasize=256;
	cv::Image<int> accum(thetasize, rhosize);
	cvZero( &accum );

	// convert ridges into lines in rho, theta space
	for(int i=0; i<r.size(); i++){
		float rho;
		CvPoint p;
		float theta;
		r.get(i, &p, &theta);
		if(theta>M_PI) theta-=M_PI;

		theta-=M_PI/2;
		rho = p.x*cos(theta) + p.y*sin(theta);
		if(rho < 0){
			rho = -rho;
			theta -= M_PI;
		}
		while(theta<0) theta+=2*M_PI;
			
		printf("x=%d y=%d rho=%f, %d theta=%f,%f rhosize=%d\n", p.x, p.y, rho, cvFloor(rho), theta, M_PI*2*cvFloor(theta/M_PI*.5*thetasize)/thetasize, rhosize);
		assert(rho<rhosize);
		assert(rho>=0);
		//LineDetector::drawLine(gimage, gscale*rho, theta, cvScalar(255,255,0)); 
		//LineDetector::drawLine(gimage, gscale*cvFloor(rho), M_PI*2*cvFloor(theta/M_PI*.5*thetasize)/thetasize, cvScalar(0,255,0)); 
		//cvShowImage("win", gimage);
		//cvWaitKey();
		accum[ cvFloor(rho) ][ cvFloor(theta/M_PI*.5*thetasize) ]++;
	}
	
	for(int j=0; j<accum.height; j++){
		for(int i=0; i<accum.width; i++){
			if(accum[j][i]>=2){
				float rho = j;
				float theta = i*M_PI*2/accum.width;
				lines.push_back( cvPoint2D32f(rho, theta) );
			}
		}
	}
}

int main(int argc, char ** argv){
	cvRedirectError( cvSegReport);
	cv::Image<uchar, 3> im, dst;
	cv::Image<float, 1> imf;
	cv::Image<uchar, 3> im3;
	std::vector< CvPoint2D32f > clusters;
	RidgeDetector ridgedetector;
	cvNamedWindow("win", 0);
	cvNamedWindow("res", 0);
	for(int i=1; i<argc; i++){
		assert( im.open( argv[i] ) );
		im3 = im;
		gimage = &im3;
		im3.show("win");
		ridgedetector.detect_multiscale( &im );
		draw( &im3, ridgedetector, 1, CV_RGB(255, 0 ,0) );
		clusters.clear();
		cluster(ridgedetector, L2(cvPoint2D32f(0,0), cvPoint2D32f(im.width, im.height)), clusters);
		draw_clusters( &im3, clusters, 1, CV_RGB(0,0,255) );
		im3.show("res");
		cvWaitKey();
	}
	cvSetMouseCallback("win", cb, NULL);
	while(1){
		im3 = im;
		float besttheta=0;
		float bestv=0;
		for(int i=0; i<8; i++){
			float theta = i*M_PI*2/8;
			float v = ridgedetector.measure( gx, gy, theta );
			if(v > bestv){
				bestv = v;
				besttheta = theta;
			}
		}
		ridgedetector.draw(&im3, CV_RGB(255,0,0));
		cvLine(&im3, cvPoint(gx,gy), cvPoint(gx+10*cos(besttheta), gy+10*sin(besttheta)), CV_RGB(0,255,0), 1);
		im3.show("win");
	}
}

