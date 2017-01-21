#ifndef __BULLSEYE_DETECTOR_HH
#define __BULLSEYE_DETECTOR_HH
#include <cv/blob_detector.hh>


//! class for detecting specific colored blobs in an image

class BullsEyeDetector : public BlobDetector{
protected:
	cv::Image<uchar> m_threshImage;
	int m_center;
	int m_areaThreshold;
	
public:
    //! Constructor
    /**
     * Construct a blob detector that detects all colors between lower and upper
     * @param lower the low end of the color range allowed by this blob detector
     */
    BullsEyeDetector(): BlobDetector(), m_center(-1), m_areaThreshold(4) {
	}

    
	bool isInside(CvConnectedComp * inblob, CvConnectedComp * outblob){
		if(inblob->rect.x >= outblob->rect.x &&
		   inblob->rect.y >= outblob->rect.y &&
		   inblob->rect.x + inblob->rect.width <= outblob->rect.x + outblob->rect.width &&
		   inblob->rect.y + inblob->rect.height <= outblob->rect.y + outblob->rect.height){
		   	return true;
		}
		return false;
	}
    /// detect blobs
    /**
     *  @param img binary input image
     *  @param in_place ignored by this class
     */
    virtual void detect(const IplImage *img){
		// first threshold img
		m_threshImage.reshape(img);
		cvAdaptiveThreshold(img, &m_threshImage, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY_INV, 21);
		
		//matlab_compat::figure(2);
		//matlab_compat::imshow(&m_threshImage);

		BlobDetector::detect(&m_threshImage, false);
		
		// reject any blobs that aren't square or have a small area
		for(int i=m_blobs->total-1; i>=0; i--){
			CvConnectedComp * blob = this->getBlob(i);
			// allow n pixel tollerance
			double ratio = blob->rect.width / (double) blob->rect.height;
			if(ratio < 1) ratio = 1/ratio;
			if(ratio > 2 || blob->area<=m_areaThreshold){
				cvSeqRemove(m_blobs, i);
			}
		}
		
		// return the blob that encloses 2 other blobs
		for(int i=0; i<m_blobs->total-1; i++){
			CvConnectedComp * blob = this->getBlob(i);
			int ninside = 0;
			
			for(int j=i+1; j<m_blobs->total; j++){
				CvConnectedComp * cblob = this->getBlob(j);
				
				// test if I'm inside
				if(isInside(cblob, blob)){
					ninside++;
				}
			}
			if(ninside==2){
				m_center = i;
				return;
			}
		}

		// none found
		this->clear();
		m_center=-1;

	}
	CvConnectedComp * getCenter() {
		return m_center==-1 ? NULL : this->getBlob(m_center);
	}
	void draw(IplImage * im, const CvScalar & color){
		if(m_blobs->total==0) return;
		CvConnectedComp * blob = this->getCenter();
		CvPoint center = cvPoint(blob->rect.x+blob->rect.width/2, blob->rect.y+blob->rect.height/2);
		int radius = blob->rect.width/2;
		cvCircle( im, center, radius, color );
	}
};

#endif //__BULLSEYE_DETECTOR_HH
