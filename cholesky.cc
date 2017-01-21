#include <cv/matrix.hh>
#include <cv/cvext.h>

int main(int argc, char ** argv){
	Matrix<float> mat(10,10), truth(10,10), cholesky(10,10);
	for(int i=0; i<mat.rows; i++){
		for(int j=0; j<mat.cols; j++){
			if(i>=j){
				mat[i][j] = drand48();
			}
			else{
				mat[i][j] = 0;
			}
		}
	}
	truth = mat*mat.t();

	cvCholesky( &truth, &cholesky );
	std::cout<<mat<<std::endl;
	std::cout<<truth<<std::endl;
	std::cout<<cholesky<<std::endl;
	cholesky = cholesky*cholesky.t();
	std::cout<<cholesky<<std::endl;
}
