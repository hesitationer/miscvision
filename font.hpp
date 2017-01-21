#ifndef __CV_FONT_HPP
#define __CV_FONT_HPP
#include <cxcore.h>

namespace cv {

class Font : public CvFont {
public:
	Font(double scale=1.0) {
		cvInitFont(this, CV_FONT_HERSHEY_PLAIN, scale, scale, 0, 1, CV_AA);
	}
	int height(){
		CvSize sz;
		int basey;
		cvGetTextSize( "Fyg", this, &sz, &basey);
		return sz.height+basey;
	}
	void render(CvArr * arr, const char * text, const CvScalar & color){
		cvPutText(arr, text, cvPoint(0, this->height()), this, color);
	}
};

}
#endif //__CV_FONT_HPP
