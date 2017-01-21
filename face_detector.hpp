#ifndef __FACE_DETECTOR_HH
#define __FACE_DETECTOR_HH

#include <cv.h>
#include <cv/haar_detector.hpp>

/// class to detect human faces in an image - uses adaboost/haar-wavelet implementation in opencv
class FaceDetector : public HaarDetector{
public:
	FaceDetector(const char * cascade_fname) :
		HaarDetector(cascade_fname)
	{
	}
	
	FaceDetector() :
		HaarDetector("../data/haarcascades/haarcascade_frontalface_alt.xml")
	{
	}

	/** return the number of faces found in the last filtration */
	int getNumFaces(){
		return this->size();
	}
	/** return the ith face found in the last filtration */
	CvRect * getFace(int i){
		return (CvRect*)cvGetSeqElem( m_objects, i);
	}
};

#endif //__FACE_DETECTOR_HH
