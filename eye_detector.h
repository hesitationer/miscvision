#ifndef __EYE_DETECTOR_H
#define __EYE_DETECTOR_H

#include <cv/image.h>
#include <cv/gaussian.h>
#include <cv/dog.h>
#include <cv/face_detector.h>

class EyeDetector {
protected:
	CvRect m_right;
	CvRect m_left;
	int m_width;  // size of subsampled face image
	cv::Image<float> m_lpdf;
	cv::Image<float> m_rpdf;
	cv::Image<float> m_imf;
	cv::Image<float> m_dog;
	cv::Image<float> m_temp;

public:
	EyeDetector():
		m_width(96)
			//m_right(cvRect(0,0,0,0)),
			//m_left(cvRect(0,0,0,0))
	{
	}

	CvRect getRightEye() { return m_right; }
	CvRect getLeftEye() { return m_left; }

	/// assume given intensity image is of head region only
	void detect(const cv::Image<uchar> & intensity, const CvRect & roi, bool flipped=false){
		if(m_lpdf.reshape(m_width,m_width)){
			calcEyePdf(m_lpdf, true, flipped);
		}
		if(m_rpdf.reshape(m_width,m_width)){
			calcEyePdf(m_rpdf, false, flipped);
		}


		m_imf.reshape(m_width,m_width);
		m_imf.convert(intensity.scale(m_width, m_width));
		m_dog.reshape(m_imf);
		m_temp.reshape(m_imf);

		cvDog(&m_imf, &m_dog, &m_temp, m_width/3+1, m_width/3+1);
		
		// left eye
		CvPoint minp,maxp;
		double minv,maxv;
		m_temp = m_dog * m_lpdf;
		cvMinMaxLoc(&m_temp, &minv, &maxv, &minp, &maxp);
		maxp.x = roi.x + intensity.width*(maxp.x + .5 - m_width/8 )/m_width;
		maxp.y = roi.y + intensity.height*(maxp.y + .5 - m_width/8 )/m_width;
		m_left = cvRect(maxp.x, maxp.y, (m_width/4)*intensity.width/m_width, (m_width/4)*intensity.height/m_width);
		
		// right eye
		m_temp = m_dog * m_rpdf;
		cvMinMaxLoc(&m_temp, &minv, &maxv, &minp, &maxp);
		maxp.x = roi.x + intensity.width*(maxp.x + .5 - m_width/8 )/m_width;
		maxp.y = roi.y + intensity.height*(maxp.y + .5 - m_width/8 )/m_width;
		m_right = cvRect(maxp.x, maxp.y, (m_width/4)*intensity.width/m_width, (m_width/4)*intensity.height/m_width);
	
	}

	void calcEyePdf(cv::Image<float> &map, bool left, bool flipped){
		Gaussian x1,y1;
		int w=map.width;
		int h=map.height;

		if(left){
			x1.setMean(w*(1-1/1.4));
		}
		else{
			x1.setMean(w*(1/1.4));
		}
		x1.setSig(0.03*w);
		y1.setMean(h*(flipped ? 1/1.6 : 1-1/1.6));
		y1.setSig(0.02*h);

		float scale = x1.pdf(x1.getMean())*y1.pdf(y1.getMean());
		for(int j=0;j<h; j++){
			for(int i=0; i<w; i++){
				map[j][i] = x1.pdf(i)*y1.pdf(j)/scale;
			}
		}
	}

	void draw(IplImage * im, CvScalar color){
		cvCircle(im, cvPoint(m_left.x+m_left.width/2, m_left.y+m_left.height/2), m_left.width/2, color);
		cvCircle(im, cvPoint(m_right.x+m_right.width/2, m_right.y+m_right.height/2), m_right.width/2, color);
		
		//cvRectangle(im, cvPoint(m_left.x, m_left.y), cvPoint(m_left.x+m_left.width, m_left.y+m_left.height), color);
		//cvRectangle(im, cvPoint(m_right.x, m_right.y), cvPoint(m_right.x+m_right.width, m_right.y+m_right.height), color);
	}
};

#endif //__EYE_DETECTOR_H
