#ifndef __CORNER_DETECTOR_HPP
#define __CORNER_DETECTOR_HPP

#include <vector>
#include <cv/fixed_priority_queue.hpp>
#include <cv.h>
#include <cv/image.h>
#include <cv/cvext.hpp>

class CornerDetector {
	cv::Image<uchar, 1> m_ismax;
	cv::Image<float, 1> m_response;
public:
	typedef enum {
		HARRIS,
		EIGEN,
		SOBEL
	} Method;

	void detect( CvArr * arr, std::vector< CvPoint > & features, int num=20, int wsize=3, int method=EIGEN, CvArr * prob=NULL){
		CvSize sz = cvGetSize(arr);
		m_ismax.realloc( sz.width, sz.height );
		m_response.realloc( sz.width, sz.height );
		features.clear();

		switch (method) {
		case HARRIS:
			cvCornerHarris( arr, &m_response, wsize);
			break;
		case EIGEN:
			cvCornerMinEigenVal( arr, &m_response, wsize );
			break;
		case SOBEL:
			cvPreCornerDetect( arr, &m_response, wsize );
			break;
		default:
			printf("Unrecognized corner detection method given");
			return;
		}
		fixed_priority_queue< CvPoint > queue( num );
		cvIsLocalMax( m_response, m_ismax );
		if(prob){
			//m_response = *cv::Image<float,1>::safe_cast( prob );
			cvMul( &m_response, prob, &m_response );
		}
		for(int j=0; j<m_ismax.height; j++){
			for(int i=0; i<m_ismax.width; i++){
				if(m_ismax[j][i]!=0 && m_response[j][i]>0){
						queue.push( cvPoint(i,j), m_response[j][i] );
				}
			}
		}
		while(!queue.empty()){
			features.push_back( queue.top().data );
			//printf("%d %d : %f\n", queue.top().data.x, queue.top().data.y, queue.top().weight);
			queue.pop();
		}
		std::reverse(features.begin(), features.end());
	}
	const cv::Image<float,1> response(){
		return m_response;
	}
	const cv::Image<uchar,1> maxima(){
		return m_ismax;
	}
};

#endif // __CORNER_DETECTOR_HPP
