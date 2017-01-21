#ifndef __HISTOGRAM_HH
#define __HISTOGRAM_HH

#include <cxtypes.h>
#include <cxcore.h>
#include <cv.h>
#include "matrixnd.hpp"

namespace cv {
class Histogram : public CvHistogram{
public:
    Histogram(const int & dims, const int & dim_size ){
        int * dim_sizes = new int[dims];
        for(int i=0; i<dims; i++){
            dim_sizes[i] = dim_size;
        }
        init(dims, dim_sizes);
		delete dim_sizes;
    }
	Histogram(const int & dims, const int * dim_sizes){
		init(dims, dim_sizes);
	}
	Histogram( const cv::Histogram & hist ){
		int dims;
		int size[CV_MAX_DIM];
		CvHistogram * myptr = this;
		dims = cvGetDims( hist.bins, size);

		this->init( dims, size );
		cvCopyHist( &hist, &(myptr) );
	}
    void init(const int & dims, const int * dim_sizes ){
		this->thresh2 = 0;
		this->type = CV_HIST_MAGIC_VAL;
		this->bins = cvInitMatNDHeader( &this->mat, dims, dim_sizes, CV_32F);
		cvCreateData( this->bins );
    }
	void calc(CvArr ** images, bool accumulate=false, const CvArr * mask=NULL){
		if(!accumulate){double min, max;
			int dims = this->mat.dims;
			float * range = new float[2*dims];
			float ** ranges = new float*[dims];
			for(int i=0; i<dims; i++){
				cvMinMaxLoc(images[i], &min, &max);
				range[i*2] = 0;
				range[i*2+1] = max;
				ranges[i] = range+(i*2);
			}
			cvSetHistBinRanges( this, ranges);
			delete [] range;
			delete [] ranges;
		}
		cvCalcHist((IplImage **) images, this, accumulate, mask);
	}
	void calc(CvArr * image, bool accumulate=false, const CvArr * mask=NULL){
		CvArr * planes[]={image, 0, 0, 0};
		int type = cvGetElemType( image );
		int nch = CV_MAT_CN( type );
		int depth = CV_MAT_DEPTH( type );
		if(nch>1){
			for(int i=0; i<nch; i++){
				CvSize sz = cvGetSize(image);
				planes[i] = cvCreateMat( sz.height, sz.width, depth );
			}
			cvSplit(image, planes[0], planes[1], planes[2], planes[3] );
		}

		this->calc(planes, accumulate, mask);
		
		if(nch>1){
			for(int i=0; i<nch; i++){
				cvReleaseMat( (CvMat **) (planes+i) );
			}
		}
	}
    double compare(const Histogram & h){
        return cvCompareHist(this, &h, CV_COMP_CORREL);
    }
	void normalize( double factor=1 ){
		cvNormalizeHist( this, factor );
	}
	void draw(IplImage * im, CvScalar color=cvScalarAll(255), int line=CV_FILLED){
		// get dimensions from ROI
		CvRect roi = cvGetImageROI(im);
		float max_value;

		cvGetMinMaxHistValue( this, 0, &max_value, 0, 0 );
		if(max_value==0) max_value = 1;

		int nbins = this->mat.dim[0].size;
		float hscale = (float)roi.width/nbins;
		
		if(line==CV_FILLED){
			cvZero(im);
		}
		for(int h = 0; h < nbins; h++ )
		{
			float bin_val = cvQueryHistValue_1D( this, h);
			float intensity = bin_val/max_value;
			//CvScalar color = cvScalarAll(h*255/nbins); 
			cvRectangle( im, cvPointFrom32f( cvPoint2D32f(h*hscale, roi.height-1 )),
					cvPointFrom32f( cvPoint2D32f((h+1)*hscale - 1, roi.height - intensity*roi.height - 1)),
					color, line );
		}
	}
	void stretch(IplImage * im, IplImage * res){
		CvScalar avg, std;
		CvMat mbins;
		cvGetMat(this->bins, &mbins, NULL, 1);
		cvAvgSdv(&mbins, &avg, &std);
	
		int lidx,uidx;
		for(lidx=0; lidx<this->mat.dim[0].size; lidx++){
			if(cvQueryHistValue_1D(this, lidx) > avg.val[0] ) break;
		}
		for(uidx=this->mat.dim[0].size-1; uidx>=0; uidx--){
			if(cvQueryHistValue_1D(this, uidx) > avg.val[0] ) break;
		}

		double uthresh,lthresh;
		uthresh = uidx*this->thresh[0][1]/(this->mat.dim[0].size); 
		lthresh = lidx*this->thresh[0][1]/(this->mat.dim[0].size); 

		cvThreshold(im, res, uthresh, 255, CV_THRESH_TRUNC);
		
		// the lower threshold is slightly tricky -- 
		// first shift all values so that lthresh is at zero
		cvSubS(res, cvScalarAll(lthresh), res);
		// now threshold at 0
		cvThreshold(res,res, 0, 255, CV_THRESH_TOZERO);

		// shift back
		cvAddS(res, cvScalarAll(lthresh), res);
	}
	bool save(const char * fname){
		cvSave( fname, this );
		return true;
	}

	bool open(const char * fname){
		CvHistogram * hist = (CvHistogram *)cvLoad(fname);
		assert(CV_IS_UNIFORM_HIST(hist));

		cvDecRefData(&(this->mat));

		memcpy(this->thresh, hist->thresh, sizeof(this->thresh));
		memcpy(&(this->mat), &(hist->mat), sizeof(CvMatND));

		cvFree(&hist);
		return true;
	}

};

inline std::ostream & operator <<(std::ostream & out, const Histogram & h) {
	switch(h.mat.dims){
		case 1:
			return out<<*(cv::MatrixND<float, 1> *)h.bins;
		case 2:
			return out<<*(cv::MatrixND<float, 2> *)h.bins;
		case 3:
			return out<<*(cv::MatrixND<float, 3> *)h.bins;
		default:
			break;
	}
	return out;
}

} // namespace cv

#endif //__HISTOGRAM_HH
