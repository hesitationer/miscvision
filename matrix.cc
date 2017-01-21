#include <stdio.h>
#include <cv/cvext.h>
#include <cv/matrix.hpp>

cv::Matrix<float, 1> * newmat(){
	cv::Matrix<float, 1> mat(10,10);
	cv::Matrix<float, 1> * nmat = new cv::Matrix<float,1>(mat, -1, 0);
	return nmat;
}
void mul( const cv::Matrix<float> & A, const cv::Matrix<float> & B, cv::Matrix<float> & C){
	cvMatMul(&A, &B, &C);
}

int main(int argc, char ** argv){
	cvRedirectError( cvSegReport );
	cv::Matrix<float, 1> A(10000,600), B(600,1), C(10000,1);
	cv::Matrix<float, 1> * mat = new cv::Matrix<float,1> (10,10);
	mul(A, B, C);
	for(int i=0; i<mat->rows; i++){
		for(int j=0; j<mat->cols; j++){
			(*mat)[i][j] = i*mat->cols+j;
		}
	}
	for(int i=0; i<mat->rows; i++){
		cv::Matrix<float, 1> row = mat->row(i);
		for(int j=0; j<row.cols; j++){
			printf("%f ", row[0][j]);
		}
		printf("\n");
	}

	cv::Matrix<float, 1> * m1, * m2;
	m1 = new cv::Matrix<float, 1>(*mat);
	m2 = new cv::Matrix<float, 1>(*mat, cvRect(0,0, 5, 5));

	delete mat;

    for(int i=0; i<m1->rows; i++){
        cv::Matrix<float, 1> row = m1->row(i);
        for(int j=0; j<row.cols; j++){
            printf("%f ", row[0][j]);
        }
        printf("\n");
    }
	
	delete m1;

    for(int i=0; i<m2->rows; i++){
        cv::Matrix<float, 1> row = m2->row(i);
        for(int j=0; j<row.cols; j++){
            printf("%f ", row[0][j]);
        }
        printf("\n");
    }
	m2->save("A.yml");

	delete m2;

	m1 = newmat();
	printf("%d x %d : %f\n", m1->rows, m1->cols, (*m1)[0][0]);

	m1->load("A.yml");
    for(int i=0; i<m1->rows; i++){
        cv::Matrix<float, 1> row = m1->row(i);
        for(int j=0; j<row.cols; j++){
            printf("%f ", row[0][j]);
        }
        printf("\n");
    }

	delete m1;
}
