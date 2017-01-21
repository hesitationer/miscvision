#ifndef __CV_ELLIPSE_DETECTOR_HPP
#define __CV_ELLIPSE_DETECTOR_HPP

#include <cv/image.h>
#include <cv/seq.hpp>
#include <cv/contour_detector.hpp>

void cvEllipseBox32f(IplImage * im, CvBox2D *box, CvScalar color,
		             int thickness CV_DEFAULT(1),
		             int line_type CV_DEFAULT(8), int shift CV_DEFAULT(0) ){
	CvPoint center;
	CvSize size;
	// Convert ellipse data from float to integer representation.
	center.x = cvRound(box->center.x);
	center.y = cvRound(box->center.y);
	size.width = cvRound(box->size.width*0.5);
	size.height = cvRound(box->size.height*0.5);

	//fprintf(stdout,"angle: %f\n", box->angle);
	// Draw ellipse.
	cvEllipse(im, center, size,
			-box->angle, 0, 360,
			color, thickness, line_type, shift);
}

int cmp_box_area(const void* _a, const void* _b, void* _userdata){
	CvBox2D32f *a=(CvBox2D32f*)_a,*b=(CvBox2D32f*)_b;
	return -(a->size.width*a->size.height - b->size.width*b->size.height);
}

class EllipseDetector : public ContourDetector {
protected:
	cv::Seq<CvBox2D32f> m_ellipses;
	CvMemStorage * m_storage;
public:
	EllipseDetector():ContourDetector(){
		m_storage = cvCreateMemStorage(0);
		m_ellipses.init(m_storage);
	}
	void detect(IplImage * im, int thresh){
		CvBox2D32f box;
		CvPoint* PointArray;
		CvPoint2D32f* PointArray2D32f;
		ContourDetector::detect(im, thresh);

		// This approximate each contour by an ellipse.
		for(CvSeq* cont=m_contours;cont;cont = cont->h_next)
		{   
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

			// Get contour point set.
			cvCvtSeqToArray(cont, PointArray, CV_WHOLE_SEQ);

			// Convert CvPoint set to CvBox2D32f set.
			for(i=0; i<count; i++)
			{
				PointArray2D32f[i].x = (float)PointArray[i].x;
				PointArray2D32f[i].y = (float)PointArray[i].y;
			}

			// Fits ellipse to current contour.
			cvFitEllipse(PointArray2D32f, count, &box);

			m_ellipses.push(box);
			//cvSeqPush(m_ellipses, &box);
		}
		cvSeqSort( &m_ellipses, cmp_box_area );
	}
	CvBox2D32f getEllipse(int i){
		return m_ellipses[i];
		//(CvBox2D32f *) cvGetSeqElem(m_ellipses, i);
	}
	CvBox2D32f & operator[](int i){
		return m_ellipses[i];
	}
	int getNumEllipses(){
		return m_ellipses.size();//->total;
	}
	int size(){
		return m_ellipses.size();
	}
	void clear(){
		ContourDetector::clear();
		//cvClearSeq(m_ellipses);
		m_ellipses.clear();
	}
};

#endif //__CV_ELLIPSE_DETECTOR_HPP
