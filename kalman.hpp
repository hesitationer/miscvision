#ifndef __CV_KALMAN_HPP
#define __CV_KALMAN_HPP

#include <cxcore.h>
#include <cv.h>

namespace cv {

template <typename model_t, typename measurement_t>
class Kalman : public CvKalman {
public:
	Kalman() {
		int ndp, nmp;
		ndp = sizeof(model_t)/sizeof(float);
		nmp = sizeof(measurement_t)/sizeof(float);
		CvKalman * k = cvCreateKalman( ndp, nmp);
		memcpy(this, k, sizeof(CvKalman));
		cvFree( (void **) &k);
	}
	~Kalman() {
		// this is a little silly, but unfortunately neccessary
		CvKalman * k = (CvKalman*) cvAlloc( sizeof(CvKalman) );
		memcpy(k, this, sizeof(CvKalman));
		cvReleaseKalman( &k );
	}
	model_t * predict(){
		return (model_t *) cvKalmanPredict( this );
	}

	model_t * correct( const measurement_t * measurement){
		CvMat mat;
		cvInitMatHeader( &mat, this->MP, 1, CV_32F, (float *) measurement);
		return (model_t *) cvKalmanCorrect( this, &mat )->data.ptr; 
	}
	model_t * getPredictedState(){
		return (model_t *) this->state_pre->data.ptr;
	}
	model_t * getPreState(){
		return (model_t *) this->state_pre->data.ptr;
	}
	model_t * getPostState(){
		return (model_t *) this->state_post->data.ptr;
	}
};

} // namespace cv

#endif // __CV_KALMAN_HPP
