#ifndef __LINE_FILTER_HH
#define __LINE_FILTER_HH

#include <math.h>
#include <iostream>
using namespace std;
#include <cxcore.h>
#include <cv.h>

/// opencv based hough line detector
class LineDetector {
protected:
	CvMemStorage * m_lineStorage;
	CvSeq * m_lineSeq;
	int m_width;
	int m_height;
public:
	LineDetector(){
		m_lineSeq=0;
		m_lineStorage=0;
	}
	
    void detect(CvArr * im, bool in_place = false){
		CvSize sz;
		sz = cvGetSize(im);
		if(!m_lineStorage){
			m_lineStorage = cvCreateMemStorage(0);
		}
		if(m_lineSeq){
			cvClearMemStorage(m_lineStorage);
		}
		m_width = sz.width;
		m_height = sz.height;
		/*m_lineSeq = cvHoughLines2(edge, m_lineStorage, 
									CV_HOUGH_PROBABILISTIC, 1, CV_PI/180, 100, 
									10, 10 );
*/
       // m_lineSeq = cvHoughLines2( edge, m_lineStorage, CV_HOUGH_STANDARD, 1, CV_PI/2, 1, 10, 10 );
//        m_lineSeq = cvHoughLines2( m_image, m_lineStorage, CV_HOUGH_STANDARD, 1, 2*CV_PI, 1, 10, 10 );
        m_lineSeq = cvHoughLines2( im, m_lineStorage, CV_HOUGH_STANDARD, 4, CV_PI/16, 1, .5, 0.25);
	}
	int size(){
		return m_lineSeq->total;	
	}
	CvPoint * operator[] (int i){
		return (CvPoint *)cvGetSeqElem(m_lineSeq,i);
	}
    void getLine(int i, CvPoint2D32f * pt1, CvPoint2D32f * pt2){
		float * line = (float *) cvGetSeqElem(m_lineSeq, i);
		assert( LineDetector::fromHoughLine(line[0], line[1], cvSize(m_width, m_height), pt1, pt2) );
    }
	static void toHoughLine( CvPoint2D32f p1, CvPoint2D32f p2, float *rho, float *theta){
		*theta = atan2( p2.x-p1.x, p2.y-p1.y );
		*rho = p2.x * cos(*theta) + p2.y * sin(*theta);
	}
	static bool fromHoughLine( float rho, float theta, CvSize imsize, CvPoint2D32f * p1, CvPoint2D32f * p2 ){
        double a = cos(theta), b = sin(theta);
        cout<<rho<<","<<theta<<endl;
        if( fabs(b) < 0.001 )
        {
                p1->x = p2->x = rho;
                p1->y = 0;
                p2->y = imsize.height-1;//color_dst->height;
        }
        else if( fabs(a) < 0.001 )
        {
                p1->y = p2->y = rho;
                p1->x = 0;
                p2->x = imsize.width-1;//color_dst->width;
        }
        else
        {
            CvPoint2D32f p[4];
            p[0].x=0;
            p[0].y = LineDetector::solveY(rho, theta, p[0].x);

            p[2].y=0;
            p[2].x = LineDetector::solveX(rho, theta, p[2].y);

            p[3].y=imsize.height-1;
            p[3].x = LineDetector::solveX(rho, theta, p[3].y);

            p[1].x=imsize.width-1;
            p[1].y = LineDetector::solveY(rho, theta, p[1].x);
            bool one = true;
            bool two = true;
            //pick out valid points
            for(int i=0; i<2; i++){
				cout<<p[i].x<<","<<p[i].y<<endl;
                if(p[i].y>=0 && p[i].y<imsize.height){
                    if(one){
                        *p1=p[i];
                        one = false;
                    }
                    else{
                        *p2=p[i];
                        two = false;
                    }
                }
				cout<<p[i+2].x<<","<<p[i+2].y<<endl;
                if(p[i+2].x >=0 && p[i+2].x < imsize.width){
                    if(one){
                        *p1=p[i+2];
                        one = false;
                    }
                    else{
                        *p2=p[i+2];
                        two = false;
                    }
                }
            }

            if(one || two){
				return false;
			}

            //sort by x
            if(p1->x > p2->x){
                CvPoint2D32f t = *p1;
                *p1 = *p2;
                *p2 = t;
            }
        }
		return true;			
	}
    static float solveX(float rho, float theta, float y){
        //rho = xcos theta + y sin theta
        return (rho - y*sin(theta))/cos(theta);
    }
    static float solveY(float rho, float theta, float x){
        return (rho - x*cos(theta))/sin(theta);
    }
	void draw(IplImage * im, CvScalar color){
		for(int i=0; i<this->size() && i<10; i++){
			float * line = (float *) cvGetSeqElem(m_lineSeq, i);
			this->drawLine(im, line[0], line[1], color);
		}
	}
    static void drawLine(IplImage * im, float rho, float theta, CvScalar color, int thickness=1){
        CvPoint2D32f p1, p2;
		//if(LineDetector::fromHoughLine(rho, theta, cvGetSize(im), &p1, &p2)){
		//	printf("%f %f -> [%f,%f] [%f,%f]\n", rho, theta, p1.x, p1.y, p2.x, p2.y);
		double a = cos(theta), b = sin(theta);
            double x0 = a*rho, y0 = b*rho;
            p1.x = (x0 + 1000*(-b));
            p1.y = (y0 + 1000*(a));
            p2.x = (x0 - 1000*(-b));
            p2.y = (y0 - 1000*(a));
	        cvLine(im, cvPointFrom32f(p1), cvPointFrom32f(p2), color, thickness);
		//}
    }
	
};	
#endif //__FILTERED_IMAGE_HH

