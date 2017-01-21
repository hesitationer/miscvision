#include <cv/image.h>
#include <cv/matrix.hh>
#include <cv/cvext.h>
#include <cv/line_detector.hpp>
#include <cv/image_pyramid.hh>
#include <cv/ridge_detector.hpp>

class GaussianImage : public cv::Image<float> {
public:
    /** downsampes image by a factor of 2 */
    static GaussianImage & dscale2(const GaussianImage & a, GaussianImage & b) {
        int w=a.width,h=a.height;
        b.realloc(w/2,h/2);
        cvPyrDown(&a, &b);
        return b;
    }
	GaussianImage & operator=(const cv::Image<float> & im){
		cv::Image<float>::operator=(im);
		return (*this);
	}
	
};
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
IplImage * gimage;
float gscale=1.0;

void extract( cv::Image<uchar,3> & im, RidgeDetector & r, float scale, std::vector< CvRect > & rects){
	cv::Image<uchar, 3> zim;
	char key;
	for(int i=0; i<r.size(); i++){
		CvPoint p;
		float theta;
		r.get(i, &p, &theta);
		zim = im;
		p.x *= scale;
		p.y *= scale;
		//cvRotateZoomImage(&im, &zim, cvPointTo32f( p ), theta, 4*scale/zim.width);
		CvRect r = cvRect( p.x - 2*scale, p.y - 2*scale, 4*scale, 4*scale );
		cvRectangle( &zim, cvPoint(r.x, r.y), cvPoint(r.x+r.width, r.y+r.height), CV_RGB(255,128, 0), 1);
		zim.show("win");
		key = (char)cvWaitKey();
		switch(key){
		case 'y':
			rects.push_back( r );
			break;
		case 'q':
			return;
		case 'k':
			return;
		default:
			continue;
		}
	}
}

int main(int argc, char ** argv){
	cvRedirectError( cvSegReport);
	cv::Image<uchar, 3> im, dst;
	cv::Image<float, 1> imf;
	cv::Image<uchar, 3> im3;
	cv::ImagePyramid< GaussianImage > pyr;
	std::vector< CvPoint2D32f > clusters;
	std::vector< CvRect > rects;
	RidgeDetector ridgedetector;
	
	FILE * positiveF, * negativeF;
	positiveF = fopen(argv[1], "w");
	negativeF = fopen(argv[2], "w");

	assert(positiveF);
	assert(negativeF);

	cvNamedWindow("win", 0);
	cvNamedWindow("res", 0);
	for(int i=3; i<argc; i++){
		assert( im.open( argv[i] ) );
		imf = im;
		im3 = im;
		gimage = &im3;
		im3.show("win");
		pyr.calculate( imf );
		rects.clear();
		for(int j=0; j<pyr.size(); j++){
			gscale = pow(2.0,j);
			ridgedetector.detect( &(pyr[j]) );
			draw( &im3, ridgedetector, pow(2.0, j), CV_RGB(255, j*255/pyr.size(),0) );
			im3.show("res");
			if((char)cvWaitKey()=='q') break;
			im3 = im;
			extract( im3, ridgedetector, gscale, rects );
		}
		if( rects.size() > 0 ){
			fprintf(positiveF, "%s %d", argv[i], rects.size());
			for(int j=0; j<rects.size(); j++){
				CvRect r = rects[j];
				fprintf(positiveF, " %d %d %d %d", r.x, r.y, r.width, r.height);
			}
			fprintf(positiveF, "\n");
		}
		else{
			fprintf(negativeF, "%s\n", argv[i]);
		}
		if((char)cvWaitKey()=='q') break;
	}
	fclose( negativeF );
	fclose( positiveF );
}

