#include <cv/cvext.h>
#include <matrix.hh>

int main(int argc, char ** argv){
	cvRedirectError(cvSegReport);
	for(int i=1; i<argc; i++){
		Matrix<float> mat;
		mat.open(argv[i]);
		std::cout<<mat<<std::endl;
	}
}
