#include <cxcore.h>
#include <types.hpp>
template <typename pixel_t, int nch=1>
struct SparseMat {
	CvSparseMat * hdr;

	SparseMat():hdr(NULL){
	}

	SparseMat( int dims, int * sizes ){
		hdr = cvCreateSparseMat(dims, sizes, CV_MAKETYPE(CvDepthType<pixel_t>::depth(), nch));
	}
};

template <typename pixel_t, int nch=1>
struct MatND {
	CvMatND * hdr;
	MatND( int dims, int * sizes ){
		hdr = cvCreateMatND( dims, sizes, CV_MAKETYPE(CvDepthType<pixel_t>::depth(), nch));
	}
	pixel_t operator[](size_t idx){
		
	}
};

int main(int argc, char ** argv){
	int sizes[] = {15,15,15};
	SparseMat<float> mat( 3, sizes ); 
}

