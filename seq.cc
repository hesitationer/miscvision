#include <iostream>
#include <cxcore.h>
#include <algorithm>
#include <seq.hpp>

#define NPOINTS 20
struct CvPoint2D32f_CMP {
	bool operator () ( const CvPoint2D32f & a, const CvPoint2D32f & b ){
		return a.x < b.x;
	}
};

std::ostream & operator<< ( CvPoint2D32f p, std::ostream & out ) {
	out << "(" <<p.x<<","<<p.y<<")";
	return out;
}

template <class T> struct print : public std::unary_function<T, void>
{
	print(std::ostream& out) : os(out), count(0) {}

	void operator() (T x) { os << ' '; ++count; }
	std::ostream& os;
	int count;
};

int main(int argc, char ** argv){
	CvMemStorage * mem = cvCreateMemStorage(0);
	cv::Seq<CvPoint2D32f> seq(mem);

	for(int i=0; i<NPOINTS; i++){
		seq.push_back( cvPoint2D32f( drand48(), drand48() ) );
	}
	
	print< CvPoint2D32f > P = for_each( seq.begin(), seq.end(), print<CvPoint2D32f>(std::cout) );
	std::cout<< std::endl << P.count << " objects printed."<<std::endl;

	std::random_shuffle( seq.begin(), seq.end() );

	std::sort( seq.begin(), seq.end(), CvPoint2D32f_CMP() );

	for(int i=0; i<NPOINTS; i++){
		printf("%f %f\n", seq[i].x, seq[i].y);
	}
}
