#ifndef CV_EXT_HPP
#define CV_EXT_HPP

#include <cxcore.h>
#include <cv.h>
#include <cv/image.h>
#include <cv/matrix.hh>
#include <cv/cvarr_adapter.hpp>
#include <cv/fixed_priority_queue.hpp>


/////////////  Supplemental function for C++ data structures
template<typename T>
void cvMinMaxTopN( const cv::Image<T> & im, int num, float & min, float & max, CvPoint & minp, CvPoint & maxp, bool invert=false){
    fixed_priority_queue< CvPoint > queue(num);
	float scale = (invert ? -1.0 : 1.0);
	CvRect roi = im.getImageROI();
	roi.width+=roi.x;
	roi.height+=roi.y;
    for(int j=roi.y; j<roi.height; j++){
        for(int i=roi.x; i<roi.width; i++){
            //printf("%f  ?=  ", im[j][i]);
            //print_queue(queue);
            queue.push(cvPoint(i,j), scale*im[j][i]);
        }
    }
    //printf("----------\n");
    min = queue.top().weight;
	minp = queue.top().data;
    queue.pop();
    while(!queue.empty()){
        max = queue.top().weight;
		maxp = queue.top().data;
        //printf("%f\n", min);
        queue.pop();
    }
	min *= scale;
	max *= scale;

	if(min > max){
		float tmp = min;
		min = max;
		max = tmp;
	}

}
template<typename T, int nch>
void cvMinMaxTopN( const cv::Image<T,nch> & im, int num, float & min, float & max, CvPoint & minp, CvPoint & maxp, bool invert=false){ 
	printf("cvTopN not implemented for channels > 1\n");
}

template <typename T>
bool cvIsLocalMin(const cv::Image<T, 1> &im, int x, int y, int w, int h, T thresh, bool testcenter){
	int y1=MAX(0, y-h/2);
	int y2=MIN(im.height, y+h/2+1);
	int x1=MAX(0, x-w/2);
	int x2=MIN(im.width, x+h/2+1);

	for(int j=y1; j<y2; ++j){
		const T* src = im[j];
		if(j==y  && !testcenter){
			for(int i=x1; i<x; ++i){
				if(src[i]<=thresh) return false;
			}
			for(int i=x+1; i<x2; ++i){
				if(src[i]<=thresh) return false;
			}
		}
		else{
			for(int i=x1; i<x2; ++i){
				if(src[i]<=thresh) return false;
			}
		}
	}
	return true;
}

template <typename T>
bool cvIsLocalMax(const cv::Image<T, 1> &im, int x, int y, int w, int h, T thresh, bool testcenter){
	int y1=MAX(0, y-h/2);
	int y2=MIN(im.height, y+h/2+1);
	int x1=MAX(0, x-w/2);
	int x2=MIN(im.width, x+h/2+1);

	for(int j=y1; j<y2; ++j){
		const T* src = im[j];
		if(j==y && !testcenter){
			for(int i=x1; i<x; ++i){
				if(src[i]>=thresh) return false;
			}
			for(int i=x+1; i<x2; ++i){
				if(src[i]>=thresh) return false;
			}
		}
		else{
			for(int i=x1; i<x2; ++i){
				if(src[i]>=thresh) return false;
			}
		}
	}
	return true;
}

/*///////////////////////////////////////////////////////////////////////////////////////////////////////
// @param im
// @param cmp
// @param res
// @param test_center
// @param size
// @param mask
//////////////////////////////////////////////////////////////////////////////////////////////////////M*/
template <typename T, int nch>
void cvIsLocalMax2(const cv::Image<T, nch> &im, const cv::Image<T, nch> & cmp, cv::Image<uchar, 1> & res, bool test_center=false, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	printf("cvIsLocalMax2 not implemented for multi-channel matrices\n");
	assert(0);
}

template <typename T, int nch>
void cvIsLocalMin2(const cv::Image<T, nch> &im, const cv::Image<T, nch> & cmp, cv::Image<uchar, 1> & res, bool test_center=false, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	printf("cvIsLocalMin2 not implemented for multi-channel matrices\n");
	assert(0);
}

template <typename T>
void cvIsLocalMax2(const cv::Image<T, 1> &im, const cv::Image<T, 1> & cmp, cv::Image<uchar, 1> & res, bool test_center=false, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	if(!mask){
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (cvIsLocalMax(cmp, i, j, size, size, src[i], test_center) ? 255 : 0);
			}
		}
	}
	else{
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			const uchar *mptr = (*mask)[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (mptr[i]!=0 && cvIsLocalMax(cmp, i, j, size, size, src[i], test_center)? 255 : 0);
			}                                               
		} 
	}
}
template <typename T>
void cvIsLocalMin2(const cv::Image<T, 1> &im, const cv::Image<T, 1> & cmp, cv::Image<uchar, 1> & res, bool test_center=false, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	if(!mask){
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (cvIsLocalMin(cmp, i, j, size, size, src[i], test_center) ? 255 : 0);
			}
		}
	}
	else{
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			const uchar *mptr = (*mask)[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (mptr[i]!=0 && cvIsLocalMin(cmp, i, j, size, size, src[i], test_center) ? 255 : 0 );
			}                                               
		} 
	}
}

template <typename T>
void cvIsLocalMin(const cv::Image<T, 1> &im, cv::Image<uchar, 1> & res, int size=3, const cv::Image<uchar, 1> * mask=NULL ){
	int width = MIN(size, im.width);
	int height = MIN(size, im.height);
	if(!mask){
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = cvIsLocalMin(im, i, j, width, height, src[i], false) ? 255 : 0;
			}
		}
	}
	else{
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			const uchar * mptr = (*mask)[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (mptr[i]!=0 && cvIsLocalMin(im, i, j, width, height, src[i], false)) ? 255 : 0;
			}
		}
	}
}

template <typename T>
void cvIsLocalMax(const cv::Image<T, 1> &im, cv::Image<uchar, 1> & res, int size=3, const cv::Image<uchar, 1> * mask=NULL ){
	int width = MIN(size, im.width);
	int height = MIN(size, im.height);
	if(!mask){
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			for(int i=0; i<im.width; i++){
				dst[i] = cvIsLocalMax(im, i, j, width, height, src[i], false) ? 255 : 0;
			}
		}
	}
	else{
		for(int j=0; j<im.height; j++){
			const T * src = im[j];
			uchar * dst = res[j];
			const uchar * mptr = (*mask)[j];
			for(int i=0; i<im.width; i++){
				dst[i] = (mptr[i]!=0 && cvIsLocalMax(im, i, j, width, height, src[i], false)) ? 255 : 0;
			}
		}
	}
}


template <typename T, int nch>
void cvIsLocalMax(const cv::Image<T, nch> &im, cv::Image<uchar, 1> & res, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	printf("cvIsLocalMax not implemented for multi-channel matrices\n");
	assert(0);
}

template <typename T, int nch>
void cvIsLocalMin(const cv::Image<T, nch> &im, cv::Image<uchar, 1> & res, int size=3, const cv::Image<uchar, 1> * mask=NULL){
	printf("cvIsLocalMin not implemented for multi-channel matrices\n");
	assert(0);
}


template <typename data_t>
void adaptive_threshold( const cv::Image<data_t> &src, cv::Image<data_t> &dest, int size){
	cvSmooth( &src, &dest, CV_GAUSSIAN, size, size );
	int i,j;
	int rows = src.height;
	int cols = src.width;
    for( i = 0; i < rows; i++ )
    {
        const data_t * s = src[i];
        const data_t * m = dest[i];
        data_t * d = dest[i];

        for( j = 0; j < cols; j++ ){
            d[j] = s[j] - m[j];
		}
    }
}

template <typename data_t>
void cvBayerSplit(const cv::Image<data_t, 1> & im,
                  cv::Image<data_t, 1> & q1,
                  cv::Image<data_t, 1> & q2,
                  cv::Image<data_t, 1> & q3,
                  cv::Image<data_t, 1> & q4){
    q1.realloc(im.width/2, im.height/2);
    q2.realloc(im.width/2, im.height/2);
    q3.realloc(im.width/2, im.height/2);
    q4.realloc(im.width/2, im.height/2);

    for(int i=0; i<im.height; i++){
        const data_t * p = im[i];
        data_t * pq1, *pq2;
        if(i%2==0){
            pq1=q1[i/2];
            pq2=q2[i/2];
        }
        else{
            pq1=q3[i/2];
            pq2=q4[i/2];
        }
        for(int j=0; j<im.width; j+=2){
            pq1[j/2] = p[j];
			pq2[j/2] = p[j+1];
        }
    }
}

// solve for matrix L s.t. A = L*L^T  where L is lower triangular and L^T is upper triangular
// formula is as follows
//
// L(i,j) = 1/L(i,j) * ( A(i,j) - SUM(k=0...j-1, L(i,k)*L(j,k)) )  for i > j
// L(i,i) = sqrt( A(i,i) - SUM(k=0...i-1, L(i,k)*L(i,k)) ) for i==j
template <typename T>
int cvCholesky( const Matrix<T, 1> & src, Matrix<T, 1> & dest ){
	
	assert(src.rows==src.cols);
	assert(src.rows==dest.rows);
	assert(src.cols==dest.cols);

	for(int i=0; i<src.rows; i++){
		const T * srow = src[i];
		T * drow_i = dest[i];
		for(int j=0; j<src.cols; j++){
			T * drow_j = dest[j];
			if(i<j){
				drow_i[j]=0;
			}
			else{
				drow_i[j] = srow[j];

				for(int k=0; k<j; k++){
					drow_i[j] -= (drow_i[k]*drow_j[k]);
				}
				
				if(i==j){
					drow_i[j] = (T) sqrt( (double) (drow_i[j]) );
				}
				else {
					drow_i[j] /= drow_j[j];
				}
			}
		}
	}
	return 0;
}

template <typename T, int nch>
int cvCholesky( const Matrix<T, nch> & src, Matrix<T, nch> & dest ){
	printf("cvCholesky not implemented for multi-channel matrices\n");
	assert(0);
	return -1;
}

template <typename T>
double cvCrossCorrelation(const Matrix<T, 1> & src1, const Matrix<T, 1> & src2){
	// xcorr = SUM[ (x(i) - m_x) * (y(i)-m_y) ] / sqrt( SUM[ (x(i)-m_x)^2 ] * SUM[ (y(i)-m_y)^2 ] )
	
	// m_x, m_y 
	CvScalar mean1 = cvAvg( &src1 );
	CvScalar mean2 = cvAvg( &src2 );
	
	double sum1 = 0, sum2=0, sum3=0;
	int w=src1.cols;
	int h=src1.rows;
	
	for( int j=0; j<h; j++ ){
		const T * row1 = src1[j];
		const T * row2 = src2[j];
		for(int i=0; i<w; i++){
			sum1 += (row1[i]-mean1.val[0])*(row1[i]-mean1.val[0]);
			sum2 += (row2[i]-mean2.val[0])*(row2[i]-mean2.val[0]);
			sum3 += (row1[i]-mean1.val[0])*(row2[i]-mean2.val[0]);
		}
	}
	
	//  sqrt( SUM[ (x(i)-m_x)^2 ] * SUM[ (y(i)-m_y)^2 ] )
	double bottom = sqrt( sum1 * sum2 ); // cvNorm( tmp1, NULL, CV_L2 ) * cvNorm( tmp2, NULL, CV_L2 );

	return sum3/bottom;
}

#endif //__CV_EXT_HPP
