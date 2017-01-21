#ifndef __STEERED_RIDGE_HH
#define __STEERED_RIDGE_HH

#include <cv/image.h>
#include <cv/cvext.h>
#include <cv/fixed_priority_queue.hpp>

class RidgeDetector{
protected:
	cv::Image<float, 1> m_response;
	cv::Image<float, 1> m_theta;
	cv::Image<float, 1> m_dxx;
	cv::Image<float, 1> m_dyy;
	cv::Image<float, 1> m_dxy;
	cv::Image<float, 1> m_y;
	cv::Image<float, 1> m_x;
	cv::Image<uchar, 1> m_ismax;
	std::vector<CvPoint> m_maxima;
public:
	RidgeDetector(){
	}
	void detect_multiscale( IplImage * im ){
		m_response.realloc( *im );
		m_dxx.realloc( *im );
		m_dxy.realloc( *im );
		m_dyy.realloc( *im );
		m_theta.realloc( *im );
		m_ismax.realloc( *im );
		m_x.realloc( *im );
		m_y.realloc( *im );

		int minpyrsize = 16;
		// ridge function is now multiscale
		// I(x,y,theta,s) = SUM abs( sin*sin*Ixx(x,y,s) + cos*cos*Iyy(x,y,s) - 2*sin*cos*Ixy(x,y,s) ) -
		//                      abs( cos*cos*Ixx(x,y,s) + sin*sin*Iyy(x,y,s) + 2*sin*cos*Ixy(x,y,s) )
		// dI/dtheta = 4*cos(x)*sin(x)*a - 4*sin(x)*cos(x)*b - 4*c*(cos(x)^2 - sin(x)^2)
		// where 
		// a = SUM_s Ixx(x,y,s)
		// b = SUM_s Iyy(x,y,s)
		// c = SUM_s Ixy(x,y,s)
		//
		// so again
		// x = atan( 2*c/(a-b))*.5
		CvSize sz = cvGetSize(im);
		int levels = cvFloor(log( MIN(sz.width,sz.height)/minpyrsize ));
		cv::Image<float> ** pyr = (cv::Image<float> **) cvCreateImagePyr( sz, IPL_DEPTH_32F, 1, levels);
		cv::Image<float> ** tmp = (cv::Image<float> **) cvCreateImagePyr( sz, IPL_DEPTH_32F, 1, levels );
		*(pyr[0]) = *im;
		
		cvSobel(pyr[0], &m_dxx, 2, 0, 5);
		cvSobel(pyr[0], &m_dyy, 0, 2, 5);
		cvSobel(pyr[0], &m_dxy, 1, 1, 5);
		for(int i=1; i<levels; i++){
			cvPyrDown( pyr[i-1], pyr[i]);
			cvSobel(pyr[i], tmp[i], 2, 0, 5);
			cvResize( tmp[i], tmp[0], CV_INTER_NN);
			cvAdd(&m_dxx, tmp[0], &m_dxx);
			cvSobel(pyr[i], tmp[i], 0, 2, 5);
			cvResize( tmp[i], tmp[0], CV_INTER_NN );
			cvAdd(&m_dyy, tmp[0], &m_dyy);
			cvSobel(pyr[i], tmp[i], 1, 1, 5);
			cvResize( tmp[i], tmp[0], CV_INTER_NN );
			cvAdd(&m_dxy, tmp[0], &m_dxy);
		}
		cvReleaseImagePyr( (IplImage***) &pyr );
		cvReleaseImagePyr( (IplImage***) &tmp );

		this->findmaxima();
	}

	void detect(IplImage * im, bool in_place=false){
		m_response.realloc( *im );
		m_dxx.realloc( *im );
		m_dxy.realloc( *im );
		m_dyy.realloc( *im );
		m_theta.realloc( *im );
		m_ismax.realloc( *im );
		m_x.realloc( *im );
		m_y.realloc( *im );

		// convert to floating point
		m_response = *im;

		//cvSmooth(&m_response, &m_response, CV_GAUSSIAN, 23, 23);
		// compute second derivatives
		cvSobel(&m_response, &m_dxx, 2, 0, 5);
		cvSobel(&m_response, &m_dyy, 0, 2, 5);
		cvSobel(&m_response, &m_dxy, 1, 1, 5);

	 	this->findmaxima();   
	}

	void findmaxima(){

		// Now find orientation for each pixel
		// orientation is maximum theta of function 
		// I = abs( sin*sin*Ixx + cos*cos*Iyy - 2*sin*cos*Ixy ) -
		//     abs( cos*cos*Ixx + sin*sin*Iyy + 2*sin*cos*Ixy )
		// I = abs( sin(theta)^2*Ixx + cos(theta)^2*Iyy - 2*sin(theta)*cos(theta)*Ixy ) -
		//     abs( cos(theta)^2*Ixx + sin(theta)^2*Iyy + 2*sin(theta)*cos(theat)*Ixy )
		// a = Ixx, b= Iyy, c=Ixy
		// solve for theta D(I(theta), theta) = 0
		// I is not continuous, so split into 4 functions and differentiate
		// 1. D(I) = 2*cos(x)*sin(x)*a - 2*sin(x)*cos(x)*b - 2*c(cos(x)^2 - sin(x)^2) 
		//          -( -2*sin(x)*cos(x)*a + 2*cos(x)*sin(x)*b + 2*c*( cos(x)^2 - sin(x)^2 ))
		//         = 4*cos(x)*sin(x)*a - 4*sin(x)*cos(x)*b - 4*c( cos(x)^2 - sin(x)^2 )
		// 2. D(I) = 2*cos(x)*sin(x)*a - 2*sin(x)*cos(x)*b - 2*c(cos(x)^2 - sin(x)^2) +
		//           -2*sin(x)*cos(x)*a +2*cos(x)*sin(x)*b + 2*c*( cos(x)^2 - sin(x)^2 )
		//         = 0
		// 3. D(I) = -( 2*cos(x)*sin(x)*a - 2*sin(x)*cos(x)*b - 2*c(cos(x)^2 - sin(x)^2) )
		//           -( -2*sin(x)*cos(x)*a +2*cos(x)*sin(x)*b + 2*c*( cos(x)^2 - sin(x)^2 ))
		//         = 0
		// 4. D(I) = -( 2*cos(x)*sin(x)*a - 2*sin(x)*cos(x)*b - 2*c(cos(x)^2 - sin(x)^2) ) +
		//            -2*sin(x)*cos(x)*a +2*cos(x)*sin(x)*b + 2*c*( cos(x)^2 - sin(x)^2 )
		//         = -4*cos(x)*sin(x)*a +4*sin(x)*cos(x)*b + 4*c*(cos(x)^2 - sin(x)^2)
		//
		// So, D(I) = 4*cos(x)*sin(x)*a - 4*sin(x)*cos(x)*b - 4*c*(cos(x)^2 - sin(x)^2)
		// where meaningfully differentiable
		
		// simplifying ... 
		// using cos(x)^2 - sin(x)^2 = cos(2x)
		// 4*sin(x)*cos(x)*(a-b) - 4*c( cos( 2x ) ) = 0
		// using sin(2x)=sin(x)cos(x) and 
		// 2*sin(2x)*(a-b) - 4*c*cos( 2x ) = 0
		// 2*c/(a-b) = sin(2x)/cos(2x)
		// 2*c/(a-b) = tan(2x)
		// x = atan( 2*c/(a-b) )*.5
		
		cvScale( &m_dxy, &m_y, 2.0, 0 );
		cvSub( &m_dxx, &m_dyy, &m_x );
	//	m_x.imagesc("win"); cvWaitKey(-1);
	//	m_y.imagesc("win"); cvWaitKey(-1);
		cvCartToPolar( &m_x, &m_y, NULL, &m_theta);
		cvScale( &m_theta, &m_theta, 0.5 );
	//	m_theta.imagesc("win"); cvWaitKey(-1);

		//cvZero( &m_theta );
		// Now calculate response at maximum theta
		
		cv::Image<float> cos_theta;
		cv::Image<float> sin_theta;
		cos_theta.realloc( m_x );
		sin_theta.realloc( m_x );
		cvPolarToCart( NULL, &m_theta, &cos_theta, &sin_theta); 

		cv::Image<float> sincos_theta = cos_theta*sin_theta;
		cos_theta = cos_theta*cos_theta;
		sin_theta = sin_theta*sin_theta;

		cv::Image<float> r1 = sincos_theta*m_dxy*-2;
		cv::Image<float> r2 = r1*-1;
		
		// sin^2*dxx + cos^2*dyy - 2*sin*cos*dxy
		m_x = sin_theta*m_dxx;
		m_y = cos_theta*m_dyy;
		r1 = r1 + m_x;
		r1 = r1 + m_y;

		cvAbs(&r1, &r1);

		// cos^2*dxx + sin^2*dyy + 2*sin*cos*dxy
		m_x = cos_theta*m_dxx;
		m_y = sin_theta*m_dyy;
		r2 = r2 + m_x;
		r2 = r2 + m_y;

		cvAbs(&r2, &r2);
	//	m_response.imagesc("win"); cvWaitKey(-1);
	//	m_x.imagesc("win"); cvWaitKey(-1);
	//	m_y.imagesc("win"); cvWaitKey(-1);
		m_response = r1 - r2;

	//	m_response.imagesc("win"); cvWaitKey(-1);

		cvIsLocalMax(m_response, m_ismax);
		// Now find maxima
		fixed_priority_queue< CvPoint > queue(10);
		m_maxima.clear();
		for(int j=0; j<m_response.height; j++){
			for(int i=0; i<m_response.width; i++){
				//print_queue(queue);
				if(m_ismax[j][i]>0 && m_response[j][i]>0){
					queue.push(cvPoint(i,j), m_response[j][i]);
				}
			}
		}
		queue.pop();
		while(!queue.empty()){
			CvPoint p = queue.top().data;
			//printf("%d,%d = %f, %f \n", p.x, p.y, m_theta[p.y][p.x], m_response[p.y][p.x]);
			m_maxima.push_back( queue.top().data );
        	queue.pop();
    	}
    }
	void get(int i, CvPoint * p, float * theta){
		*p = m_maxima[i];
		*theta = m_theta[p->y][p->x];
	}
	void draw(CvArr * arr, CvScalar color){
		float rad = 10;
		for(int i=0; i<m_maxima.size(); i++){
			CvPoint p = m_maxima[i];
			CvPoint p1, p2;
			float theta = m_theta[p.y][p.x];
			p1 = cvPoint(p.x - rad*cos(theta), p.y - rad*sin(theta));
			p2 = cvPoint(p.x + rad*cos(theta), p.y + rad*sin(theta));
			cvCircle(arr, p, 1, color, 1, CV_AA);
			cvLine(arr, p1, p2, color, 1, CV_AA);
		}

	}
	int size() const {
		return m_maxima.size();
	}

    double measure(int x, int y, double theta){
        CvScalar dxx = cvGet2D(&m_dxx, y, x);
        CvScalar dyy = cvGet2D(&m_dyy, y, x);
        CvScalar dxy = cvGet2D(&m_dxy, y, x);
        return fabs(sin(theta)*sin(theta)*dxx.val[0] + 
               cos(theta)*cos(theta)*dyy.val[0] - 
               2*sin(theta)*cos(theta)*dxy.val[0]) -
			  fabs( sin(theta)*sin(theta)*dyy.val[0] +
					cos(theta)*cos(theta)*dxx.val[0] +
					2*sin(theta)*cos(theta)*dxy.val[0] );
		// orientation is maximum theta of function 
		// I = sin(theta)^2*Ixx + cos(theta)^2*Iyy - 2*sin(theta)*cos(theta)*Ixy
		// 2*cos(x)*sin(x)*a - 2*sin(x)*cos(x)*b + 2*c(cos(x)^2 - sin(x)^2) = 0
		// 2*sin(x)*cos(x)*(a-b) + 2*c( cos( 2x ) ) = 0
		// sin(2x)*(a-b) + 2*c*cos( 2x ) = 0
		// 2*c/(b-a) = sin(2x)/cos(2x)
		// 2*c/(b-a) = tan(2x)
		// x = atan( 2*c/(b-a) )*.5
    }
};

#endif //__EDGE_DETECT_HH
