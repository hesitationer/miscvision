#include <cv/image.h>
#include <cv.h>

namespace cv {

CvPoint match_template( const CvArr * intensity, const CvArr * tmplt, CvRect window, int method=CV_TM_CCORR){ 
	// 1. Figure out search window
	CvMat roi;
	CvSize isize = cvGetSize(intensity);
	CvSize tsize = cvGetSize(tmplt);
	
	//std::cout<<"match_template("<<window.x<<","<<window.y<<","<<window.width<<","<<window.height<<")"<<std::endl;
	
	// fix window to be within image dimensions
	if(window.x<0) window.x=0;
	if(window.y<0) window.y=0;
	if(window.width+window.x>=isize.width) window.width=isize.width-window.x-1;
	if(window.height+window.y>=isize.height) window.height=isize.height-window.y-1;

	// ensure template fits in window
	if(window.width < tsize.width) return cvPoint(-1,-1);//assert(0);
	if(window.height < tsize.height) return cvPoint(-1,-1);//assert(0);

	//std::cout<<"mt fixed args:("<<window.x<<","<<window.y<<","<<window.width<<","<<window.height<<")"<<std::endl;


	// 2. Cross correlate template and intensity in given window 
	IplImage * result = cvCreateImage(cvSize(window.width-tsize.width+1,
	            window.height-tsize.height+1), IPL_DEPTH_32F, 1);
	//cv::Image<float> result(window.width-tsize.width+1, 
	//		window.height-tsize.height+1);
	cvGetSubRect(intensity, &roi, window);
	cvMatchTemplate( &roi, tmplt, result, method);

	// 3. Find the best match
	CvPoint minp,maxp;
	double minv,maxv;
	cvMinMaxLoc( result, &minv, &maxv, &minp, &maxp);

	if(method==CV_TM_SQDIFF_NORMED || method==CV_TM_SQDIFF){
		return cvPoint( window.x + minp.x ,  window.y + minp.y );
	}
	return cvPoint( window.x + maxp.x ,  window.y + maxp.y );

}
CvPoint match_template( const CvArr * intensity, const CvArr * tmplt, const CvPoint & pos, CvSize radius, int method = CV_TM_CCORR){ 
	CvSize sz = cvGetSize( tmplt );
	return match_template( intensity, tmplt, cvRect(pos.x-radius.width, pos.y-radius.height, sz.width+radius.width*2, sz.height+radius.height*2), method);
}


}
