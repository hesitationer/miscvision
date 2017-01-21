#include <cxcore.h>
#include <cv/match_template.h>

namespace cv {

class StereoImage {
public:
	IplImage * leftImage;
	IplImage * rightImage;
	float cz;
	int maxDisparity;
	
	StereoImage(IplImage * left, IplImage * right, float a_cz, int a_maxDisparity):
		leftImage( left ),
		rightImage( right ),
		cz( a_cz ),
		maxDisparity( a_maxDisparity ) 
	{
	}
	static int calcDisparity(const IplImage * left, const IplImage * right, int x, int y, int maxDisparity=-1, int radius=10){
		// 1. Define template
		CvMat tmplt;
		CvRect roi = cvRect(x-radius, y-radius, radius*2, radius*2);
		if(roi.x<0) roi.x=0;
		if(roi.y<0) roi.y=0;
		if(roi.x+roi.width>left->width) roi.width = left->width - roi.x;
		if(roi.y+roi.height>left->height) roi.height = left->height - roi.y;

		if(maxDisparity==-1){
			maxDisparity = roi.x;
		}

		cvGetSubRect( left, &tmplt, roi);
		
		// 2. find matching point in right image within m_maxDisparity
		CvPoint p = cv::match_template(right, &tmplt, cvRect(roi.x-maxDisparity, roi.y, roi.width+maxDisparity, roi.height));
		if(p.x==-1) return 0;

		return  roi.x - p.x;
	}
	int calcDisparity(int x, int y, int radius){
		return cv::StereoImage::calcDisparity(leftImage, rightImage, x, y, maxDisparity, radius);
	}
	static float calcDepth(float disparity, float cz){
		return cz/disparity;
	}
	float calcDepth(float disparity) {
		return this->cz/disparity;
	}
};

}
