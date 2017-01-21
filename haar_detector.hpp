#ifndef __HAAR_DETECTOR_HH
#define __HAAR_DETECTOR_HH

#include <cv.h>
#include <cv/image.h>

/// class to detect human objects in an image - uses adaboost/haar-wavelet implementation in opencv
class HaarDetector {
public:
/**
 * The default parameters (scale_factor=1.1, min_neighbors=3, prune=false) 
 * are tuned for accurate yet slow face detection. For faster face detection 
 * on real video images the better settings are 
 * (scale_factor=1.2, min_neighbors=2, prune=true). 
 */
	HaarDetector(const char * fname ):
		m_cascade(NULL){
    	m_storage = cvCreateMemStorage(0);
		m_object_size = cvSize(40,40);
		m_haar_scale_factor = 1.1;
		m_neighbors = 3; 
		m_prune = false;
		m_objects = NULL;
		m_resizeImage = false;
		if( !this->open(fname) )
		{
			fprintf( stderr, "ERROR: Could not load classifier cascade '%s'\n", fname);
			return;	
		}
	}
	bool open(const char * cascade_fname){
		if(m_cascade){
			cvReleaseHaarClassifierCascade( &m_cascade );
		}
		m_cascade = (CvHaarClassifierCascade*)cvLoad( cascade_fname, 0, 0, 0 );
		return m_cascade != NULL;
	}
	
	/** return the number of objects found */
	int size(){
		return (m_objects ? m_objects->total : 0);
	}

	void clear(){
		cvClearMemStorage( m_storage );
	}

	/** bracket accessor to detected objects */
	CvRect & operator[](int i){
		return *(CvRect*)cvGetSeqElem( m_objects, i);
	}

	/** set whether to down sample input image before detection or not */
	void setDownSample(bool do_downsample, int size=0){
		m_resizeImage = do_downsample;
		if(m_resizeImage){
			m_size = size;
		}
	}
	
	bool detect(const CvArr * arr){
		IplImage header;
		IplImage * im = cvGetImage( arr, &header );
		if(im->depth==8 && im->nChannels==1){
			return this->detect(*cv::Image<uchar,1>::safe_cast( im ) );
		}
		else{
			cv::Image<uchar,1> tmp;
			tmp = *im;
			return this->detect( tmp );
		}
	}

	bool detect(const cv::Image<uchar, 1> & input){
		float scale=1.0;

		// scale image to requested size
		if(m_resizeImage && m_size < input.width && m_size < input.height){
			if( input.width > input.height ){
				m_temp.realloc( input.width*m_size/input.height, m_size );
			}
			else{
				m_temp.realloc( m_size, input.height*m_size/input.width );
			}
			cvResize(&input, &m_temp, CV_INTER_NN);
			scale = input.width/(float)m_temp.width;
		}
		else{
			m_temp = input;
		}
		
		this->clear();

		m_objects = cvHaarDetectObjects( &m_temp, m_cascade, m_storage, m_haar_scale_factor, m_neighbors, (m_prune ? CV_HAAR_DO_CANNY_PRUNING:0), m_object_size );
		
		for(int i = 0; i < (m_objects ? m_objects->total : 0); i++ ){
			CvRect * r = (CvRect*)cvGetSeqElem( m_objects, i);
			r->x      = cvRound(r->x*scale);
			r->width  = cvRound(r->width*scale);
			r->y      = cvRound(r->y*scale);
			r->height = cvRound(r->height*scale);
		}

		return m_objects && m_objects->total>0;
	}

	void draw(CvArr * im, CvScalar color=CV_RGB(0,255,0)){
		for(int i=0; i<this->size(); i++){
			CvRect r = (*this)[i];
			cvRectangle(im, cvPoint(r.x, r.y), cvPoint(r.x+r.width, r.y+r.height), color, 1);
		}
	}

	virtual ~HaarDetector(){
		if(m_storage) cvReleaseMemStorage( &m_storage );
		if(m_cascade) cvReleaseHaarClassifierCascade( &m_cascade );
	}

protected:
	CvMemStorage* m_storage;
	CvHaarClassifierCascade* m_cascade;
	CvSeq * m_objects;
	cv::Image<uchar, 1> m_temp;
	bool m_resizeImage;
	double m_haar_scale_factor;
	int m_neighbors;
	int m_size;
	bool m_prune;
	CvSize m_object_size;
};

#endif //__HAAR_DETECTOR_HH
