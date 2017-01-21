#ifndef CV_CAPTURE_HPP
#define CV_CAPTURE_HPP

#include <highgui.h>

namespace cv{
	class Capture {
		CvCapture * m_capture;
	public:
		Capture():
			m_capture(NULL)
		{
		}
		Capture( const char * uri ) :
			m_capture(NULL)
		{
			this->open(uri);	
		}
		operator CvCapture * (){ return m_capture; }
		bool open( const char * uri ){
			if(!this->open_file( uri )){
				int idx = atoi(uri);
				return this->open_camera( idx );
			}
			return true;
		}
		bool open_camera( int idx ){
			m_capture = cvCreateCameraCapture( idx );
			return(m_capture!=NULL);
		}
		bool open_file( const char * filename ){
			m_capture = cvCreateFileCapture( filename );
			return (m_capture!=NULL);
		}
		~Capture(){
			cvReleaseCapture( &m_capture );
		}

	};
}

#endif  // CV_CAPTURE_HPP
