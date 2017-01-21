#ifndef __FACE_FEATURES_HPP
#define __FACE_FEATURES_HPP

#include <cv/image.h>

class FaceFeaturesPDF {
	cv::Image<uchar, 1> ** eyePDF;
	cv::Image<uchar, 1> ** facePDF;
public:
	typedef enum {
		LEYE = 0,
		REYE,
		LMOUTH,
		RMOUTH,
		LLEYEBROW,
		RLEYEBROW,
		LREYEBROW,
		RREYEBROW,
		LHEAD,
		LLEYE,
		RLEYE,
		LREYE,
		RREYE,
		RHEAD,
		NOSE,
		LNOSTRIL,
		RNOSTRIL,
		UPPERLIP,
		LOWERLIP,
		CHIN,
		MAXFEATURENUM
	} FeatureID;
	
	bool open(const char * dir){
		char fname[256];
		bool ret=true;

		// allocate arrays
		eyePDF = new cv::Image<uchar, 1> * [ FaceFeaturesPDF::MAXFEATURENUM ];
		facePDF = new cv::Image<uchar, 1> * [ FaceFeaturesPDF::MAXFEATURENUM ];

		// open images
		for(int i=0; i< FaceFeaturesPDF::MAXFEATURENUM; i++){
			sprintf(fname, "%s/face_pdf_%02d.png", dir, i);
			facePDF[i] = new cv::Image<uchar,1>;
			//printf("Opening fname %s\n", fname);
			assert(facePDF[i]->open(fname));
			sprintf(fname, "%s/eye_pdf_%02d.png", dir, i);
			eyePDF[i] = new cv::Image<uchar,1>;
			assert(eyePDF[i]->open(fname));
		}

		return ret;
	}

	cv::Image<uchar> & getFacePDF( int i ){
		return *(facePDF[i]);
	}
	
	CvPoint2D32f pointToFaceSpace(const CvRect & r, const CvPoint2D32f & p) const{
		return cvPoint2D32f( ((p.x-r.x)*facePDF[0]->width)/r.width, 
			   			     ((p.y-r.y)*facePDF[0]->height)/r.height );
	}
	
	CvPoint2D32f pointFromFaceSpace(const CvRect &r, const CvPoint2D32f & p) const{
		return cvPoint2D32f( (p.x*r.width)/facePDF[0]->width + r.x,
				             (p.y*r.height)/facePDF[0]->height + r.y );
	}

	// given face rectangle, return maximally likely position of other features.
	CvPoint2D32f findFeature(CvRect faceRect, FeatureID idx) const{
		CvPoint min,max;
		double minv, maxv;
		assert(idx < FaceFeaturesPDF::MAXFEATURENUM && idx>=0);
		cvMinMaxLoc(facePDF[idx], &minv, &maxv, &min, &max );
		return pointFromFaceSpace(faceRect, cvPointTo32f(max));
	}
	// given face rectangle, a point, and a feature number return the value of the bin
	// this point falls in
	float pdf(const CvRect &r, const CvPoint2D32f &p, FeatureID idx) const{
		float min = (1.0/255);

		if(idx>=FaceFeaturesPDF::MAXFEATURENUM) return min;

		// transform point into face space
		CvPoint2D32f q = pointToFaceSpace(r, p);
		CvPoint s = cvPointFrom32f(q);
		
		// assume pdf is minimal outside of rectangle
		if(s.x<0 || s.x >= facePDF[idx]->width ||
		   s.y<0 || s.y >= facePDF[idx]->height ){
			return min;
		}

		return (*facePDF[idx])[s.y][s.x]*(1.0/255) + min;
	}
	void getFeatures(CvRect faceRect, CvPoint2D32f * points){
		for(int i=0; i<FaceFeaturesPDF::MAXFEATURENUM; i++){
			points[i] = findFeature(faceRect, (FeatureID) i);
		}
	}
	void draw(IplImage * im, CvRect faceRect, CvScalar color){
		for(int i=0; i<FaceFeaturesPDF::MAXFEATURENUM; i++){
			cvCircle(im, cvPointFrom32f( this->findFeature( faceRect, (FeatureID) i ) ), 1, color);
		}
	}

	// given eye positions, return maximally likely position of other features.
	void findFeature(CvPoint2D32f eyes[2], CvPoint2D32f * features){
		// calculate transformation

	}
	~FaceFeaturesPDF(){
		for(int i=0; i<FaceFeaturesPDF::MAXFEATURENUM; i++){
			delete facePDF[i];
			delete eyePDF[i];
		}
	}
};
#endif //__FACE_FEATURES_HPP
