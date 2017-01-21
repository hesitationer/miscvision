#ifndef __CV_IMAGE_PYRAMID_HH
#define __CV_IMAGE_PYRAMID_HH

#include <vector>
#include <cv/wincompat.h>
#include <cv/image.h>
#include <cv/scaled_image.h>

namespace cv {

// image_t must support dscale2, realloc 
// operator=(cv::Image)
template <typename image_t> class ImagePyramid {
protected:
	std::vector<image_t> _images;
	int _logMinSize;
	int _minSize;
	image_t _exemplar;

public:	
	ImagePyramid(int minImageSize=16, image_t exemplar=image_t()):
		_logMinSize((int)log2(minImageSize)),
		_minSize(minImageSize),
		_exemplar(exemplar)
	{
			
	}
	ImagePyramid( const ImagePyramid<image_t> & pyr ):
		_images( pyr._images ),
		_logMinSize( pyr._logMinSize ),
		_minSize( pyr._minSize ),
		_exemplar( pyr._exemplar )
	{
	}

	~ImagePyramid(){
	}
	template <typename pixel_t>
	void check(const Image<pixel_t> & im){
		for(int j=0; j<im.height; j++){
			for(int i=0; i<im.width; i++){
				assert(!isnan(im(i,j)));
				assert(!isinf(im(i,j)));
				assert(im(i,j)<=255);
				assert(im(i,j)>=0);
			}
		}
	}
	template <typename pixel_t>
	void calculate(const Image<pixel_t> & im) {
		checkSizes(im);
	
		assert(_images.size()>0);

		_images[0] = im;
		for(size_t i=1; i<_images.size(); i++){
			image_t::dscale2(_images[i-1], _images[i]);
		}
	}
	template <typename pixel_t>
	void checkSizes(const Image<pixel_t> & im){
		int w,h;
		size_t size;
		
		// calculate number of scales to consider
		w = (int) log2(im.getWidth());
		h = (int) log2(im.getHeight());
		if(w>h) size = h - _logMinSize + 1;
		else size = w - _logMinSize + 1;
		//std::cout<<"imsize="<<im.getWidth()<<"x"<<im.getHeight()<<std::endl;
		//std::cout<<"minsize = "<<_minSize<<std::endl;
		//std::cout<<"logmin = "<<_logMinSize<<std::endl;
		//std::cout<<"size = "<<size<<std::endl;

		w = im.getWidth();
		h = im.getHeight();
		
		assert(size>0);
		
		if(_images.size()!=size){
			_images.resize(size, _exemplar);
		}
		
		for(size_t i=0; i<size; i++){

			assert(w>=_minSize);
			assert(h>=_minSize);
			
			_images[i].realloc(w,h);

			w /= 2;  
			h /= 2;
		}
	}
	const image_t & operator[] (int l) const{
		return _images[l];
	}
	const image_t & getImage(int l) const{
		return _images[l];
	}
	image_t & operator[] (int l){
		return _images[l];
	}
	image_t & getImage(int l){
		return _images[l];
	}
	const int & getNumScales() const{
		return _images.size();
	}
	const int & getWidth() const {
		return getImage(0).getWidth();
	}
	const int & getHeight() const {
		return getImage(0).getHeight();
	}
	int size() { return _images.size(); }

};

} //namespace cv
#endif //__CV_IMAGE_PYRAMID_HH
