#ifndef STEREO_CAMERA_HPP
#define STEREO_CAMERA_HPP

#include <cvaux.h>  // CvCalibFilter
#include <libdc1394++/DC1394Camera.hpp>
#include <cv/image.h>

class StereoCamera {
protected:
	DC1394::Camera m_leftCamera;
	DC1394::Camera m_rightCamera;
	cv::Image<uchar, 3> m_image3;
	cv::Image<uchar, 1> m_leftImage;
	cv::Image<uchar, 1> m_rightImage;
	IplImage m_leftHeader;
	IplImage m_rightHeader;
	CvCalibFilter m_calibration;
	bool m_rectify;

public:
	StereoCamera()
	{
		m_rectify = false; 
	}
	bool loadCalibration( const char * fname ){
		m_calibration.LoadCameraParams( fname );
		m_rectify = m_calibration.IsCalibrated();
		return m_rectify;
	}
	bool load(const char * lfname, const char * rfname){
		if(!m_leftImage.open(lfname) || !m_rightImage.open(rfname)) return false;

		m_image3.realloc(m_leftImage);

		if(m_rectify && m_calibration.GetStereoParams()!=NULL){
			cvCvtColor( &m_leftImage, &m_image3, CV_BayerBG2BGR );
			IplImage * srcarr[] = { &m_image3, &m_rightImage};
			IplImage * dstarr[] = { &m_image3, &m_rightImage};
			m_calibration.Rectify( srcarr, dstarr ) ;
			cvCvtColor( &m_image3, &m_leftImage, CV_BGR2GRAY ); 
		}
		return true;
	}

	StereoCamera(u_int64_t leuid, u_int64_t reuid):
		m_leftCamera(leuid),
		m_rightCamera(reuid)
	{
		assert(m_leftCamera.setup(-1, FORMAT_SVGA_NONCOMPRESSED_1, MODE_1024x768_MONO, SPEED_400, FRAMERATE_15));
		assert(m_rightCamera.setup(-1, FORMAT_SVGA_NONCOMPRESSED_1, MODE_1024x768_MONO, SPEED_400, FRAMERATE_15));
		assert(m_leftCamera.isValid());
		assert(m_leftCamera.start());
		assert(m_rightCamera.isValid());
		assert(m_rightCamera.start());
		
		cvInitImageHeader(&m_leftHeader, cvSize(m_leftCamera.getWidth(), m_leftCamera.getHeight()), 8, 1);
		cvInitImageHeader(&m_rightHeader, cvSize(m_rightCamera.getWidth(), m_rightCamera.getHeight()), 8, 1);

		this->loadCalibration("cameras.txt");
	}
	bool lock(){
		assert( m_leftCamera.capture() );
		assert( m_rightCamera.capture() );
		m_leftHeader.imageData = m_leftHeader.imageDataOrigin = (char *) m_leftCamera.getBuffer();
		m_rightHeader.imageData = m_rightHeader.imageDataOrigin = (char *) m_rightCamera.getBuffer();

		m_image3.realloc( m_leftHeader.width, m_leftHeader.height );
		m_leftImage.realloc( m_leftHeader.width, m_leftHeader.height );
		m_rightImage.realloc( m_rightHeader.width, m_rightHeader.height );

		// convert bayer pattern to bgr
		cvCvtColor( &m_leftHeader, &m_image3, CV_BayerBG2BGR );

		if(m_rectify && m_calibration.GetStereoParams()!=NULL){
			IplImage * srcarr[] = { &m_image3, &m_rightHeader};
			IplImage * dstarr[] = { &m_image3, &m_rightImage}; 
			 m_calibration.Rectify( srcarr, dstarr ) ;
		}
		else{
			m_rightImage = m_rightHeader;
		}

		cvCvtColor( &m_image3, &m_leftImage, CV_BGR2GRAY ); 

		return true;
	}
	const cv::Image<uchar,3> & color(){
		return m_image3;
	}
	const cv::Image<uchar,1> & left(){
		return m_leftImage;
	}
	const cv::Image<uchar,1> & right(){
		return m_rightImage;
	}
	bool unlock(){
		return true;
	}
	bool load( const char * fname ){
		return m_calibration.LoadCameraParams( fname );
	}


#if 0
	void rectify(){
		Vector3<float> e1,e2,e3;
		float Rdata[9];
		CvMat Rrect = cvMat(3,3,CV_32F,Rdata);

		cv::Image<uchar,1> tmp;
		tmp.realloc(m_leftImage);
		cvUndistort( &m_leftImage, &tmp, m_leftIntrinsic, m_leftDistortion ); 
		tmp.realloc(m_rightImage);
		cvUndistort( &m_rightImage, &tmp, m_rightIntrinsic, m_rightDistortion ); 

		e1 = T/e1.magnitude();
		e2 = 1/sqrt(T[0]*T[0]+T[1]*T[1]) * Vector3( -T[1], T[0], 0 );
		e3 = e1.cross(e3);

		Vector3::xstack( e1,e2,e3, Rrect );
		
		// rotate left camera so that epipole goes to infinity along horizontal (x) axis
		// same rotation to right camera
		// rotate right camera by R
		// adjust scale of both reference frames

		// build remap arrays 
		for(int j=0; j<height; j++){
			for(int i=0; i<width; i++){
				mapx[j][i] = 0;
				mapy[j][i] = 0;
			}
		}

		cvRemap();
	}
#endif
#endif
};
