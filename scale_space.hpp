#ifndef __CV_SCALE_SPACE_HH
#define __CV_SCALE_SPACE_HH

#include <cv/image.h>
#include <vector>

/// base scale space class

namespace cv {

template <typename pixel_t, int nch> 
class ScaleSpace : public std::vector< cv::Image<pixel_t, nch> >{
public:	
	typedef typename std::vector< cv::Image<pixel_t, nch> >::iterator iterator;
	ScaleSpace( int size ) : std::vector< cv::Image<pixel_t, nch> >( size )
	{
	}
	bool realloc( int width, int height ){
		bool reallocated=false;
		for(iterator iter = this->begin(); 
				iter!=this->end(); 
				++iter){
			reallocated = iter->realloc( width, height ) || reallocated;
		}
		return reallocated;
	}

	const int & getWidth() const {
		return this->begin()->width;
	}

	const int & getHeight() const {
		return this->begin()->height;
	}

	static void dscale2(const ScaleSpace<pixel_t, nch> & src, ScaleSpace<pixel_t, nch> &dest){
		assert( src.size() == dest.size() );
		for(size_t i=0; i<src.size(); i++){
			Image<pixel_t,nch>::dscale2(src[0], dest[0]);
		}
	}
	const ScaleSpace & operator=(const Image<pixel_t, nch> & im){
		(*this->begin()) = im;
		iterator it=this->begin(); 
		for( iterator it1=it+1; it1!=this->end(); it=it1, it1=it+1 ) {
			this->scale( *it, *it1 );
		}
	}
	void scale( const Image<pixel_t, nch> & src, Image<pixel_t, nch> & dest ){
	}
};

} //namespace cv
#endif //__CV_SCALE_SPACE_HH
