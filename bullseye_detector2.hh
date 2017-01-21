#ifndef BULLSEYE_DETECTOR_HH
#define BULLSEYE_DETECTOR_HH
#include "cv.h"
#include "highgui.h"

class BullsEyeDetector {
protected:
    CvMemStorage* m_storage;
    CvSeq* m_contours;
    CvBox2D32f* box;
    CvPoint* PointArray;
    CvPoint2D32f* PointArray2D32f;
	cv::Image<uchar> m_threshImage;
public:
BullsEyeDetector(){
    // Create dynamic structure and sequence.
    m_storage = cvCreateMemStorage(0);
    m_contours = cvCreateSeq(CV_SEQ_ELTYPE_POINT, sizeof(CvSeq), sizeof(CvPoint) , m_storage);
}

void detect( IplImage * img ){
    
	// first threshold img
	m_threshImage.reshape(img);
	cvAdaptiveThreshold(img, &m_threshImage, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY_INV, 21);
    
	// Threshold the source image. This needful for cvFindContours().
    //cvThreshold( img, &m_threshImage, slider_pos, 255, CV_THRESH_BINARY );
    
    // Find all contours.
    cvFindContours( &m_threshImage, m_storage, &m_contours, sizeof(CvContour), 
                    CV_RETR_TREE, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
    
}
void draw( IplImage * im, CvScalar & color){
    // This cycle draw all contours and approximate it by ellipses.
    for(CvSeq * cont = m_contours; cont != NULL ;cont = cont->h_next)
    {   
        // Draw current contour.
        cvDrawContours(im, cont, color, color, 0, 1, 8);
	}
}
void findEllipse(){
#if 0
		int i; // Indicator of cycle.
        int count = cont->total; // This is number point in contour
        CvPoint center;
        CvSize size;
        
        // Number point must be more than or equal to 6 (for cvFitEllipse_32f).        
        if( count < 6 )
            continue;
        
        // Alloc memory for contour point set.    
        PointArray = (CvPoint*)malloc( count*sizeof(CvPoint) );
        PointArray2D32f= (CvPoint2D32f*)malloc( count*sizeof(CvPoint2D32f) );
        
        // Alloc memory for ellipse data.
        box = (CvBox2D32f*)malloc(sizeof(CvBox2D32f));
        
        // Get contour point set.
        cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);
        
        // Convert CvPoint set to CvBox2D32f set.
        for(i=0; i<count; i++)
        {
            PointArray2D32f[i].x = (float)PointArray[i].x;
            PointArray2D32f[i].y = (float)PointArray[i].y;
        }
        
        // Fits ellipse to current contour.
        cvFitEllipse(PointArray2D32f, count, box);
        
        // Draw current contour.
        cvDrawContours(image04,cont,CV_RGB(255,255,255),CV_RGB(255,255,255),0,1, 8);
        
        // Convert ellipse data from float to integer representation.
        center.x = cvRound(box->center.x);
        center.y = cvRound(box->center.y);
        size.width = cvRound(box->size.width*0.5);
        size.height = cvRound(box->size.height*0.5);
        box->angle = -box->angle;
        
        // Draw ellipse.
        cvEllipse(image04, center, size,
                  box->angle, 0, 360,
                  CV_RGB(0,0,255), 1, CV_AA, 0);
        
        // Free memory.          
        free(PointArray);
        free(PointArray2D32f);
        free(box);
		#endif
    }
    
};

#endif //BULLSEYE_DETECTOR_HH
