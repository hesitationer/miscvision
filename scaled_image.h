#ifndef __CV_SCALED_IMAGE_HH
#define __CV_SCALED_IMAGE_HH

#include <cv/image.h>

namespace cv {

template <typename pixel_t> class ScaledImage : public Image<pixel_t, 1>{
protected:
	float _xscale;
	float _yscale;
public:
	ScaledImage<pixel_t> & operator=(const Image<pixel_t, 1> & im){
		_xscale=1;
		_yscale=1;
		Image<pixel_t>::operator=(im);
		return (*this);
	}
	void setScale(float xscale, float yscale){
		_xscale = xscale;
		_yscale = yscale;
	}
	pixel_t getSubPix(float x, float y) const{
		return Image<pixel_t>::getSubPix(x*_xscale, y*_yscale);
	}
	pixel_t & operator()(const int &x, const int &y){
		return Image<pixel_t>::operator()(cvFloor(x*_xscale), cvFloor(y*_yscale));
	}
	const pixel_t & operator()(const int &x, const int &y) const{
		return Image<pixel_t>::operator()(cvFloor(x*_xscale), cvFloor(y*_yscale));
	}
	Image<pixel_t, 1> & operator=(const Image_OP<pixel_t, 1> & op){
		return Image<pixel_t, 1>::operator=(op);
	}
};

}

#endif // __CV_SCALED_IMAGE_HH
