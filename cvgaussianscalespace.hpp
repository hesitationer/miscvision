#ifndef __CV_GAUSSIAN_SCALE_SPACE_HH
#define __CV_GAUSSIAN_SCALE_SPACE_HH

#include <vector>
#include <cv/image_pyramid.hh>

struct CvScaleSpace {
	IplImage ** images;
	int nLevels;
	int nSubLevels;
	float sigma;
}

CvGaussianPyr * cvCreateScaleSpace( CvSize size, int nscales, float sigma, int minSize ){
	CvGaussianPyr * pyr = (CvGaussianPyr *) cvAlloc( sizeof(CvGaussianPyr *) );
	pyr->gaussians = (CvImagePyr **) cvAlloc( sizeof(CvImagePyr *)*scales );
	pyr->scales = scales;
	pyr->sigma = sigma;
	
	for(int i=0; i<scales; i++){
		pyr->gaussians[i] = cvCreateImagePyr( size, IPL_DEPTH_32F, 1, minSize );
	}
	return pyr;
}

IplImage * cvGetScaleSpaceImage( CvScaleSpace * sp, int l, int sl ){
	return sp->images[l*sp->nSubLevels+sl];
}

float cvCalcScaleSpaceSigma(CvScaleSpace * sp, double k){
	return sp->sigma * exp2(k/(double)sp->nLevels);
}

void cvCalcScaleSpace( IplImage * im, CvScaleSpace * sp){
	for(int i=0; i<sp->nLevels; i++){
		if(i>0){
			cvResize( cvGetScaleSpaceImage( sp, i-1, -1 ), cvGetScaleSpaceImage( sp, i, 0 ), CV_INTER_NN );
		}
		for(int j=1; j<sp->nSubLevels; j++){
			cvGaussian( cvGetScaleSpaceImage( sp, i, j-1 ), cvGetScaleSpaceImage( sp, i, j ), 0, 0, sp->sigma );
		}
	}
}

void cvCalcDoGScaleSpace( IplImage * im, CvScaleSpace* sp ){
	IplImage * temp;

	// for k levels in the scale space,
	// the freqFactor between them is 2^1/k
	//
	// Normally, an octave would have scales sigma -> 2^(k-1)/k*sigma 
	//
	// However, we are typically looking for peaks in the scale space, in
	// which case, we need to search a scale above and below 
	// 
	// scale       gauss  dog   search
	// s            Y      N       N
	// s*2^1/k      Y      Y       N
	// s*2^2/k      Y      Y       Y
	// ...          Y      Y       Y
	// s*2^(k-1)/k  Y      Y       N
	
	// Specifically, the first two scales don't get searched, and
	// the last doesn't.
	// We make up for this by adding some overlap between octaves

	// ...
	// s*2^(k-1)/k  Y      Y       Y
	// s*2          Y      Y       Y   // also in octave n+1 
	// s*2*2^1/k    Y      Y       Y   // also in octave n+1
	// s*2*2^1/k    Y      Y           // also in octave n+1

	// Now s*2^n and s*2^(n+1/k) are searched (except for 2^0 and 2^1/k)
	// and s*2^(n+(k-1)/k) is searched (except for the very end of the pyramid)
	cvCalcScaleSpaceSigma(CvScaleSpace * sp, j);
	temp = cvCloneImage( im );	
	for(int i=0; i<dog->nLevels; i++){	
		// sigma*2^i
		cvSmooth( im, cvGetScaleSpaceImage( sp, i, 0 ), 0, 0, sigma );
		for(int j=1; j<dog->nSubLevels; j++){
			// sigma*2^(i+j/k)
			cvSmooth( cvGetScaleSpaceImage( sp, i, j-1), cvGetScaleSpaceImage( sp, i, j ), 0, 0, sigma );
			cvSub( cvGetScaleSpaceImage( sp, i, j-1 ), cvGetScaleSpaceImage( sp, i, j ),
				   cvGetScaleSpaceImage( dog, i, j-1) );
		}
		cvSmooth( cvGetScaleSpaceImage( sp, i, sp->nSubLevels-1 ), temp, 0, 0, sigma );
		cvSub( cvGetScaleSpaceImage( sp, i, sp->nSubLevels-1 ), temp, 
			   cvGetScaleSpaceImage( sp, i, sp->nSubLevels-1 ) );
		im = temp;
		
		if(i<dog->nLevels-1){
			// seed the next scale
			temp = cvCloneImage( cvGetScaleSpaceImage( sp, i+1, 0) );
			cvResize( im, temp, CV_INTER_NN );
			cvReleaseImage( &im );
			im = temp;
		}
	}
	cvReleaseImage( &temp );
}

namespace cv {


template <typename pixel_t> class GaussianScaleSpace {
	std::vector < cv::Image<pixel_t> > _gaussian;
	float _freqFactor;
	int _numLevels;
	float _sigma;

public:
	GaussianScaleSpace(){
	}
	GaussianScaleSpace(int levels):
		_gaussian(levels),  
		_freqFactor(stok(levels)),
		_numLevels(levels),
		_sigma(1.6)
		// for k levels in the scale space,
		// the freqFactor between them is 2^1/k
		//
		// Normally, an octave would have scales sigma -> 2^(k-1)/k*sigma 
		//
		// scale       gauss  dog   search
		// s            Y      N       N
		// s*2^1/k      Y      Y       N
		// s*2^2/k      Y      Y       Y
		// ...          Y      Y       Y
		// s*2^(k-1)/k  Y      Y       N
		
		// Specifically, the first two scales don't get searched, and
		// the last doesn't.
		// We make up for this by adding some overlap between octaves

		// ...
		// s*2^(k-1)/k  Y      Y       Y
		// s*2          Y      Y       Y   // also in octave 1
		// s*2*2^1/k    Y      Y       Y   // also in octave 1
		// s*2*2^1/k    Y      Y           // also in octave 1

		// Now s*2^n and s*2^(n+1/k) are searched (except for 2^0 and 2^1/k)
		// and s*2^(n+(k-1)/k) is searched (except for the very end of the pyramid)
		
	{
	}
	GaussianScaleSpace(const GaussianScaleSpace & s):
		_gaussian(s._gaussian.size()),
		_freqFactor(s._freqFactor),
		_numLevels(s._numLevels),
		_sigma(s._sigma)
		{
	
	}
	void init(int levels, float freqFactor=-1){
		_gaussian.resize(levels+3); // 2 for search + 1 for difference
		_freqFactor = stok(levels);
		_numLevels = levels;
		_sigma = 1.6;
		std::cout<<_gaussian.size()<<std::endl;
	}

	bool reshape(int w, int h){
		bool reshaped = false;
		for(size_t i=0; i<_gaussian.size(); i++){
			reshaped |= _gaussian[i].reshape(w,h);
		}
		return reshaped;
	}
	static int calcGaussianWidth(float sigma){
		return (((int)(sigma*3))*2 + 1);
	}
	static double stok (int s)
	{
		return pow (2.0, 1.0 / s);
	}

	static void dscale2(const GaussianScaleSpace & src, GaussianScaleSpace &dest){
		//float sigma=1.6;
		//int width = calcGaussianWidth(sigma);

		//I'll take the LAST image in src, downsample it
		//Image<pixel_t>::dscale2(src._gaussian[src._gaussian.size()-1], dest._gaussian[0]);
		Image<pixel_t>::dscale2(src._gaussian[0], dest._gaussian[0]);

		//now compute the blurs
		dest = dest._gaussian[0];
	}
	float calcSigma(double k){
		return _sigma * exp2(k/(double)_numLevels);
	}
	const GaussianScaleSpace & operator=(const Image<pixel_t> & im){
		int width;
		float sigma;
		assert(_gaussian.size()>0);

		// calculate gaussians
		for(size_t i=1; i<_gaussian.size(); i++){
			sigma = calcSigma(i);
			width = calcGaussianWidth(sigma); // 19);
			//width = nwidth-width + 1;
			//width = calcGaussianWidth(_freqFactor); // 19);
			//std::cerr<<"Sigma = "<<sigma<<std::endl;
			//std::cerr<<"Convolving with gaussian size "<<width<<std::endl;
			//cvSmooth(&(_gaussian[i-1]), &(_gaussian[i]), CV_GAUSSIAN, width, width);
			// opencv computes my sigma for me :)
			cvSmooth(&im, &(_gaussian[i]), CV_GAUSSIAN, 0, 0, sigma);
			
		}
		// first one
		if(&im != &(_gaussian[0])){
			cvSmooth(&im, &(_gaussian[0]), CV_GAUSSIAN, 0, 0, _sigma);
		}
		else{
			//std::cerr<<_gaussian[0].width<<"x"<<_gaussian[0].height<<std::endl;
			cvSmooth(&(_gaussian[0]), &(_gaussian[0]), CV_GAUSSIAN, 0, 0, _sigma);
		}

		return (*this);
	}

	Image<pixel_t> & operator[] (const int & i){
		return _gaussian[i];
	}

	const Image<pixel_t> & operator[] (const int & i) const{
		return _gaussian[i];
	}
	int size() { return _gaussian.size(); }
	const int &getWidth() const { return _gaussian[0].width; }
	const int &getHeight() const { return _gaussian[0].height; }

};

template <typename pixel_t> class GaussianPyramid : public ImagePyramid< GaussianScaleSpace<pixel_t> > {
	typedef Image<pixel_t> image_t;
public:
	GaussianPyramid(int minImageSize=16, int nSubScales=3):
		ImagePyramid< GaussianScaleSpace<pixel_t> > (minImageSize, GaussianScaleSpace<pixel_t>(nSubScales))
	{
	}
	const image_t & getImage(int i, int j) const{
		return (*this)[i][j];
	}
	image_t & getImage(int i, int j){
		return (*this)[i][j];
	}
	int getNumScales() { return ImagePyramid<GaussianScaleSpace<pixel_t> >::_images.size(); }
	int getNumSubScales() { return ImagePyramid<GaussianScaleSpace<pixel_t> >::_images[0].size(); }
	/** calculate effective sigma for this subspace */
	float getSigma(int i, int j){
		return pow(2.0, (double)i)*ImagePyramid<GaussianScaleSpace<pixel_t> >::_images[i].calcSigma(j);
	}
};

} //namespace cv
#endif //__CV_GAUSSIAN_SCALE_SPACE_HH
