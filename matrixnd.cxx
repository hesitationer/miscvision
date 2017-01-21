#include <stdio.h>
#include <drand48.h>
#include <matrixnd.hpp>
int main(){
	int sizes[]={11,30,13,5};
	CvRNG rng = cvRNG();
	cv::MatrixND<float, 4> m(sizes);
	cvRandArr( &rng, &m, CV_RAND_UNI, cvScalar(0), cvScalar(1) );
	for(int i=0; i<sizes[0]; i++){
		for(int j=0; j<sizes[1]; j++){
			for(int k=0; k<sizes[2]; k++){
				int idx[] = {i, j, k, 0};
				float * ptr = m.slice(i,j,k);
				float * fptr = (float *) cvPtrND( &m, idx);
				if(ptr!=fptr){
					printf("ptr != fptr %p %p\n", ptr, fptr);
				}

				for(int kk=0; kk<sizes[3]; kk++){
					idx[3]=kk;
					//ptr[ kk ] = drand48();
					CvScalar r = cvGetND( &m, idx );
					//printf("%d %d %d %d\n",  i, j, k, kk);
					if(r.val[0] == ptr[ kk ]) continue;
					if(r.val[0] == m[i][j][k][kk] ) continue;

					printf("cvGet=%f ptr[]=%f\n m[i][j][k][kk]=%f\n", r.val[0], ptr[ kk ], m[i][j][k][kk]);
				}
			}
		}
	}
	std::cout<<"Tests passed successfully\n";

	//std::cout<<m<<std::endl;
}
