#ifndef __CV_GAUSSIAN_SCALE_SPACE_HH
#define __CV_GAUSSIAN_SCALE_SPACE_HH

#include <vector>
#include <cv/image_pyramid.hh>

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
