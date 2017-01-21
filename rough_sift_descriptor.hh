#ifndef __ROUGH_SIFT_DESCRIPTOR_HH
#define __ROUGH_SIFT_DESCRIPTOR_HH

#include "sift_descriptor.hh"

namespace sift {

class RoughDescriptor {
	Descriptor<4> _fine;
	Descriptor<2> _medium;
	Descriptor<1> _rough;
public:
	RoughDescriptor()
		{}
	double dist(const RoughDescriptor & s, int level=0) const{
		switch(level){
		case 0:
			return _fine.dist(s._fine);
		case 1:
			return _medium.dist(s._medium);
		case 2:
			return _rough.dist(s._rough);
		}
		return -1;
	}

	// Ix and Iy are pre rotated to orientation of point
	// feature is assumed to be located at the center of image patch
	template <typename image_t>
	void calcFromPatch(const ScalePoint &p, const image_t & Ix, const image_t & Iy) {
		_fine.calcFromPatch(p, Ix, Iy);
		this->calcFromHigher(_fine, _medium);
		this->calcFromHigher(_medium, _rough);
	}
	

	void calc(ScalePoint &p, const cv::DoGPyramid<float> & dog) {
		calc(p, dog.getIx(p.level, p.subLevel), dog.getIy(p.level, p.subLevel));
	}
	
	// Ix and Iy are standard sobel images from the pyramid, not preprocessed
	template <typename image_t>
	void calc(ScalePoint &p, const image_t & Ix, const image_t & Iy){
		_fine.calc(p,Ix,Iy);
		this->calcFromHigher(_fine, _medium);
		this->calcFromHigher(_medium, _rough);
	}
	template <int i1, int i2>
	void calcFromHigher(Descriptor<i1> & d1, Descriptor<i2> & d2){
		d2 = 0.0f;
		for(int i=0; i<d1.length(); i++){
			for(int j=0; j<d1.length(); j++){
				for(int k=0; k<8; k++){
					d2(i/2, j/2, k)+= d1(i,j,k);
				}
			}
		}
		d2.normalize();
	}
};

} // namespace sift

#endif //__SIFT_DESCRIPTOR_HH
