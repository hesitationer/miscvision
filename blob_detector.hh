#ifndef __BLOB_DETECTOR_HH
#define __BLOB_DETECTOR_HH
#include <cv.h>
#include <cxcore.h>
#include <cv/image.h>

static int cmp_cc( const void* _a, const void* _b, void* userdata ){
    CvConnectedComp * a = (CvConnectedComp *)_a;
    CvConnectedComp * b = (CvConnectedComp *)_b;
    return (int)(b->area - a->area);
}

//! class for detecting specific colored blobs in an image

class BlobDetector {
protected:
    CvSeq * m_blobs;
    CvMemStorage * m_storage;
	double m_total;
	cv::Image<uchar> m_temp;
	cv::Image<uchar> m_mask;
	
public:
    //! Constructor
    /**
     * Construct a blob detector that detects all colors between lower and upper
     * @param lower the low end of the color range allowed by this blob detector
     */
    BlobDetector() {
		this->init();
	}
	virtual ~BlobDetector() {
		cvClearMemStorage(m_storage);
	}

	void init(){
		m_storage = cvCreateMemStorage(0);
        m_blobs = cvCreateSeq(0, sizeof(CvSeq), sizeof(CvConnectedComp), m_storage);
		m_total = 0;
    }
    
    /// detect blobs
    /**
     *  @param img binary input image
     *  @param in_place ignored by this class
     */
    virtual void detect(CvArr *img, bool is_binary=false, int min_size=0, CvArr *mask=NULL){
		this->clear();

		// make a copy of input
		m_temp=cv::ImageHeader<uchar,1>(img);
		if(mask) m_mask=cv::ImageHeader<uchar,1>( mask );

		CvScalar mean, sdv;	
		for(int i=1; i<m_temp.height-1; i++){
			uchar * row = m_temp[i];
			for(int j=1; j<m_temp.width-1; j++){
            	if(row[j]==0) continue;
				// find image maximum and minimum, compute average of small number of pixels around them
            	CvConnectedComp comp;
				if(!is_binary){
					cvSetImageROI( &m_temp, cvRect(j-1, i-1, 3, 3) );
					cvAvgSdv( &m_temp, &mean, &sdv );
					if(sdv.val[0] > 16) sdv.val[0]=16;
					cvResetImageROI( &m_temp );
				}
				else{
					sdv = cvScalarAll(0);
				}
				cvFloodFill( &m_temp, cvPoint(j,i), cvScalarAll(0), sdv, sdv, &comp, 8);
				//insert blob 
				if(comp.area > min_size){
					cvSeqPush(m_blobs, &comp);
					m_total += comp.area;
				}
			}
		}
        cvSeqSort(m_blobs, cmp_cc, 0);

	}
	void clear(){
		cvClearSeq(m_blobs);
		m_total = 0;
	}
	double getTotalArea(){
		return m_total;
	}
	CvConnectedComp & operator[](int i){
		return *(CvConnectedComp *) cvGetSeqElem(m_blobs, i);
	}
	int size(){
		return m_blobs->total;
	}
	void draw(IplImage * im, const CvScalar &color, int maxsize=-1){
		for(int i=0; i<m_blobs->total; i++){
			CvConnectedComp & blob = (*this)[i];
			if(blob.area < maxsize) break;
			
			cvRectangle(im, cvPoint(blob.rect.x, blob.rect.y),
			                cvPoint(blob.rect.x + blob.rect.width, blob.rect.y+blob.rect.height),
						    color);
		}
	}
};

#endif //__BLOB_DETECTOR_HH
