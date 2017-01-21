#ifndef __FLOW_TRACKER_HH
#define __FLOW_TRACKER_HH
#include <math.h>
#include <vector>

#include "cv/image.h"

#define FLOW_POINTS 100
class PointCorrespondLK{
protected:
	std::vector<CvPoint2D32f> m_points;
	std::vector<CvPoint2D32f> m_points2;
	cv::Image<uchar,1> m_pyr[2];
	cv::Image<float,1> m_temp32f1;
	cv::Image<float,1> m_temp32f2;
	int m_levels;
	std::vector<char> m_status;
	CvSize m_window;
	CvTermCriteria m_criteria;
	bool m_havePrior;
	int m_locked;
    int m_maxPoints;
	int m_nValid;
	int m_nPoints;
public:
	PointCorrespondLK(int maxpoints=64) 
	{
		m_maxPoints = maxpoints;
		m_levels = 4;
		m_window = cvSize(3,3);
		m_nValid=0;
		m_nPoints = 0;
		m_points.resize(m_maxPoints);
		m_points2.resize(m_maxPoints);
		m_status.resize(m_maxPoints);
		this->clear();
	}
	void filter(cv::Image<uchar,1> & im, cv::Image<uchar,1> & im2){
		detect(im);
		correspond(im, im2);
	}
	
	void setPoints(CvPoint2D32f * points, int npoints){
		m_points.resize(m_maxPoints);
		for(int i=0; i<npoints; i++){
			m_points[i] = points[i];	
		}
		m_nPoints = npoints;
	}
	void getPoints(CvPoint2D32f * points, int npoints){
		for(int i=0; i<npoints; i++){
			points[i] = m_points2[i];
		}
	}	
	CvPoint2D32f & operator [] (const int & i) { 
		assert(i<m_maxPoints);
		if(m_nPoints<=i) m_nPoints = i+1;
		return m_points[i]; 
	}

	int max_size() { return m_points.size(); }
	int size() { return m_points.size(); }

    void detect(cv::Image<uchar,1> & im){
		assert(im.nChannels==1);
		m_temp32f1.realloc( im );
		m_temp32f2.realloc( im );
		m_nPoints = m_maxPoints;
		m_points.resize(m_nPoints);
		cvGoodFeaturesToTrack(&im, &m_temp32f1, &m_temp32f2,
		                      &(m_points[0]), &m_nPoints, .1, 10.0);
		std::fill(m_status.begin(), m_status.end(), 1);
	}
	
	void correspond(const cv::Image<uchar,1> & im, const cv::Image<uchar,1> & im2){
		assert(im.nChannels==1);
		assert(im2.nChannels==1);

		m_pyr[0].realloc(im);
		m_pyr[1].realloc(im);
		
		m_points.resize(m_maxPoints);
		m_points2.resize(m_maxPoints);
		m_status.resize(m_maxPoints);

		if(m_nPoints > 0){
			CvTermCriteria criteria =  { CV_TERMCRIT_ITER/*|CV_TERMCRIT_EPS*/, 100, 0.1 };
			cvCalcOpticalFlowPyrLK(&im, &im2, &(m_pyr[0]), &(m_pyr[1]), 
				&(m_points[0]), &(m_points2[0]), m_nPoints, m_window, 
				m_levels, &(m_status[0]), NULL, criteria, 0);
			m_nValid=0;
			for(int i=0; i<m_nPoints; i++){
				if(this->isValid(i)){
					m_nValid++;
				}
				else{
					m_points2[i] = m_points[i];
				}
			}
		}
	} 
    void draw(IplImage * im){
        for(int i=0; i<m_nPoints; i++){
            if(isValid(i)){
                cvCircle(im, cvPointFrom32f(m_points2[i]), 2, CV_RGB(0,255,0));
                //cvLine(im, cvPointFrom32f(m_points[i]), cvPointFrom32f(m_points2[i]), CV_RGB(0,255,0));
            }
            else{
                cvCircle(im, cvPointFrom32f(m_points[i]), 2, CV_RGB(255,0,0));
            }
        }
    }

	CvPoint2D32f getFlow(int i){
		CvPoint2D32f p = m_points2[i];
		p.x -= m_points[i].x;
		p.y -= m_points[i].y;
		return p;
	}

	double getFlowMagnitude(int i){
		return sqrt((m_points2[i].x-m_points[i].x)*(m_points2[i].x-m_points[i].x) + (m_points2[i].y-m_points[i].y)*(m_points2[i].y-m_points[i].y));
	}
	double getFlowAngle(int i){
		return atan2(m_points2[i].y-m_points[i].y, m_points2[i].x-m_points[i].x);
	}
	int getNumPoints(){
		return m_nPoints;
	}
	int getNumValid(){
		return m_nValid;
	}
	bool isValid(int i){
		return m_status[i]==1;
	}
	void clear(){
		std::fill(m_status.begin(), m_status.end(), 0);
		m_nValid = 0;
	}
	const CvPoint2D32f & getPrePoint(int i){
		return m_points[i];
	}
	const CvPoint2D32f & getPostPoint(int i){
		return m_points2[i];
	}
	void update() {
		m_points = m_points2;
	}
    
    /*PointArray<double> &getPrePoints(PointArray<double> & p){
		return this->fillPoints(p, m_points);
    }
    
    PointArray<double> & fillPostPoints(PointArray<double> &p){
		return this->fillPoints(p, m_points2);
	}
	PointArray<double> & fillPoints(PointArray<double> & p, std::vector<CvPoint2D32f> & v){
		p.resize(nValid);
        for(int i=0; i<v.size(); i++){
            if(this->isValid(i)){
                p(j)=Vector3<double>(v[i].x, v[i].y, 1.0);
                j++;
            }
        }
        return p;
    }*/
    /*void filterOutliers(float p){
        std::vector<double> range;
        for(int i=0; i<m_nPoints; i++){
            if(this->isValid(i)){
                range.push_back(this->getFlowMagnitude(i));
            }
        }
        sort(range.begin(), range.end());
        for(int i=0; i<m_nPoints; i++){
            float m = this->getFlowMagnitude(i);
            if(m < range[range.size()/4] ||
               m > range[range.size()-1-range.size()/4]){
                m_status[i]=0;
            }
        }
    }*/
};

#endif //__FLOW_TRACKER_HH
