#ifndef __EYE_TRACKER_H
#define __EYE_TRACKER_H

#include <cv/image.h>
#include <cv/gaussian.h>
#include <cv/face_detector.h>
#include <cv/eye_detector.h>

class FaceEyeDetector {
protected:
	IplImage rgbHeader;
	cv::Image<uchar> m_intensity;
	FaceDetector m_faceDetector;
	EyeDetector m_eyeDetector;
public:
	FaceEyeDetector()
		//m_right(cvRect(0,0,0,0)),
		//m_left(cvRect(0,0,0,0))
	{
	}
	
	CvRect getRightEye() { return m_eyeDetector.getRightEye(); }
	CvRect getLeftEye() { return m_eyeDetector.getLeftEye(); }
	bool haveFace() { return m_faceDetector.getNumFaces() > 0; }
	CvRect getFace() { return *m_faceDetector.getFace(0); }
	int getNumFaces() { return m_faceDetector.getNumFaces(); }
	CvRect getFace(int i){ return *m_faceDetector.getFace(i); }

	// detect head and eye position -- this assumes buffer is flipped
	bool detectHeadAndEyes(unsigned char * rgb_template, int w, int h){
		cv::Image<uchar,3> * rgb;
		IplImage rgbHeader;

		// wrap an ipl image header around this buffer
		cvInitImageHeader( &rgbHeader, cvSize(w,h), 8, 3);
		rgbHeader.imageData = (char *) rgb_template;
		rgbHeader.imageDataOrigin = (char *) rgb_template;

		// cast it has a cv::Image class
		rgb = (cv::Image<uchar,3> *) &rgbHeader;

		//convert to greyscale
		m_intensity.convert(*rgb); 
		return this->detectHeadAndEyes(m_intensity, true);
	}

	/// detect position of eyes given head position
	/// all coordinates here are lower-left, not upper left like the rest of opencv
	void detectEyesFromHead(unsigned char * rgb_template, int w, int h, int roi_x, int roi_y, int roi_w, int roi_h){
		cv::Image<uchar,3> * rgb;
		cv::Image<uchar,3> smrgb(roi_w, roi_h);
		IplImage rgbHeader;
		
		// wrap an ipl image header around this buffer
		cvInitImageHeader( &rgbHeader, cvSize(w,h), 8, 3);
		rgbHeader.imageData = (char *) rgb_template;
		rgbHeader.imageDataOrigin = (char *) rgb_template;
		
		// cast it has a cv::Image class
		rgb = (cv::Image<uchar,3> *) &rgbHeader;
		
		CvRect roi = cvRect(roi_x, roi_y, roi_w, roi_h);	
		
		// convert to greyscale 
		rgb->getSubImage(roi, smrgb);
		m_intensity.convert( smrgb );
			
		this->detectEyesFromHead(m_intensity, roi, true); 
	}
	void flip(char * buf, int size){
		int s = size/2;
		char * rbuf = buf+size-1;
		char tmp;
		for(int i=0; i<s; i++){
			tmp=*buf;
			*buf=*rbuf;
			*rbuf=*buf;
			buf++;
			rbuf--;
		}
	}
	/// detect head and eye position
	bool detectHeadAndEyes(const cv::Image<uchar> & intensity, bool flipped=false){
		m_faceDetector.detect(intensity);
		if(flipped){
			this->flip(m_intensity.imageData, m_intensity.imageSize);
		}
		if(m_faceDetector.getNumFaces()>0){
			std::cout<<"Face found"<<std::endl;
			CvRect face = *m_faceDetector.getFace(0);
			cv::Image<uchar, 1> sub = intensity.getSubImage(face);
			detectEyesFromHead(sub, face, false);
			return true;
		}
		return false;
	}
	
	/// detect eye position given subimage containing only face
	/// and the original roi where the subimage exists
	void detectEyesFromHead(const cv::Image<uchar> & intensity, const CvRect & roi, bool flipped=false){
		m_eyeDetector.detect(intensity, roi, flipped);
	}
	
	/// don't use this yet.
	void track(unsigned char * rgb_template, int w, int h, int roi_x, int roi_y, int roi_w, int roi_h){

	}
};

#endif //__EYE_TRACKER_H
