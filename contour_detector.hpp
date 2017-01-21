#ifndef __CONTOUR_DETECTOR_HPP
#define __CONTOUR_DETECTOR_HPP

#include <cxcore.h>
#include <cv.h>

class ContourDetector {
	protected:
		CvMemStorage * m_storage;
		CvSeq * m_contours;
	cv::Image<uchar,1> m_image;
public:
	ContourDetector(){
	    // Create dynamic structure and sequence.
	    m_storage = cvCreateMemStorage(0);
	    //m_contours = cvCreateSeq(CV_SEQ_ELTYPE_POINT, sizeof(CvSeq), sizeof(CvPoint) , m_storage);
		m_contours = NULL;
	}
	void detect(CvArr *im, double thresh, bool binary){

		// create a local copy
		if(CV_IS_MAT(im)){
			m_image = *(CvMat *) im;
		}
		else if(CV_IS_IMAGE(im)){
			m_image = *(IplImage *) im;
		}

		if(!binary){
			cvThreshold( im, &m_image, thresh, 255, CV_THRESH_BINARY_INV );
		}

    	// Find all contours.
		cvFindContours( &m_image, m_storage, &m_contours, sizeof(CvContour), 
				        CV_RETR_LIST, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
    
	}
	void clear(){
		if(m_contours && m_contours>0) cvClearSeq(m_contours);
	}
	const cv::Image<uchar,1> & getThreshImage(){
		return m_image;
	}
	void draw(IplImage * im, CvScalar color){
		// This cycle draw all contours and approximate it by ellipses.
		for(CvSeq *cont=m_contours;cont;cont = cont->h_next)
		{
			cvDrawContours(im,cont,color, color,0,1,8,cvPoint(0,0));
		}
	}
	CvSeq * head(){
		return m_contours;
	}

};

#endif //__CONTOUR_DETECTOR_HPP
