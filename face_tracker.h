#ifndef __EYE_TRACKER_H
#define __EYE_TRACKER_H

#include <cv/image.h>
#include <cv/gaussian.h>
#include <cv/face_detector.h>

class EyeTracker {
protected:
	IplImage m_rgbHeader;
	cv::Image<uchar> m_intensity;
	FaceDetector m_faceDetector;
	EyeDetector m_eyeDetector;

public:
	EyeTracker()
	{
	}
	
	CvRect getRightEye() { return m_eyeDetector.getRightEye(); }
	CvRect getLeftEye() { return m_eyeDetector.getLeftEye(); }

	/// detect position of eyes given head position
	void detect(unsigned char * rgb_template, int w, int h, int origin=0){
		cv::Image<uchar,3> * rgb;
		IplImage rgbHeader;

		// wrap an ipl image header around this buffer
		cvInitImageHeader( &rgbHeader, cvSize(w,h), 8, 3);
		rgbHeader.imageData = (char *) rgb_template;
		rgbHeader.imageDataOrigin = (char *) rgb_template;

		// cast it has a cv::Image class
		rgb = (cv::Image<uchar,3> *) &rgbHeader;

		rgb->origin = origin;
		
		m_intensity = (*rgb); 
		this->detect(m_intensity);
	}

	/// all coordinates here are lower-left, not upper left like the rest of opencv
	void detect(unsigned char * rgb_template, int w, int h, int roi_x, int roi_y, int roi_w, int roi_h){
		cv::Image<uchar,3> * rgb;
		cv::Image<uchar,3> smrgb(roi_w, roi_h);
		
		// wrap an ipl image header around this buffer
		cvInitImageHeader( &m_rgbHeader, cvSize(w,h), 8, 3);
		m_rgbHeader.imageData = (char *) rgb_template;
		m_rgbHeader.imageDataOrigin = (char *) rgb_template;
		
		// cast it has a cv::Image class
		rgb = (cv::Image<uchar,3> *) &m_rgbHeader;
		
		CvRect roi = cvRect(roi_x, roi_y, roi_w, roi_h);	
		
		// convert to greyscale 
		rgb->getSubImage(roi, smrgb);
		m_intensity = smrgb;
		
		this->detect(m_intensity, roi); 
		
	}
	void detect(const cv::Image<uchar> & intensity){
		m_faceDetector.detect(intensity);
		if(m_faceDetector.getNumFaces()>0){
			std::cout<<"Face found"<<std::endl;
			CvRect face = *m_faceDetector.getFace(0);
			cv::Image<uchar, 1> sub = intensity.getSubImage(face);
			detect(sub, face);
		}
	}

	void detect(const cv::Image<uchar> & intensity, const CvRect & roi){
		m_eyeDetector.detect(intensity, roi);
	}
	
	void track(unsigned char * rgb_template, int w, int h, int roi_x, int roi_y, int roi_w, int roi_h){

	}
};

#endif //__EYE_TRACKER_H
