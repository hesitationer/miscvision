#include <cv/cvext.h>
#include <matrix.hh>

int main(int argc, char ** argv){
	cvRedirectError(cvSegReport);
	int dim=25;
	Matrix<float> mf(dim,dim);
	Matrix<double> df;

	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			mf[i][j] = drand48();
		}
	}

	mf.save("mf.cvm");

	Matrix<float> mf2;
	mf2.open("mf.cvm");

	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			assert(mf[i][j]==mf2[i][j]);

		}
	}

	df.open("mf.cvm");
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			assert(mf[i][j]==df[i][j]);

		}
	}

	Matrix<float> subm = mf.asshape(dim*dim);
	subm.save("mf.cvm");

	mf2.open("mf.cvm");
	for(int i=0; i<subm.rows; i++){
		for(int j=0; j<subm.cols; j++){
			assert(subm[i][j]==mf2[i][j]);
		}
	}


	printf("%d %d\n", mf.rows, mf.cols);
	subm = mf.rowrange(1, dim/5).asshape(1);
	printf("%d %d\n", subm.rows, subm.cols);
	assert(subm.rows==1);
	assert(subm.cols==dim*(dim/5-1));
	int k=0;
	for(int i=1; i<dim/5; i++){
		for(int j=0; j<dim; j++){
			assert(subm[0][k]==mf[i][j]);
			k++;
		}
	}
	subm.save("mf2.cvm");
	mf2.open("mf2.cvm");
	assert(mf2.rows==1);
	assert(mf2.cols==dim*(dim/5-1));
	for(int i=0; i<subm.rows; i++){
		for(int j=0; j<subm.cols; j++){
			assert(subm[i][j]==mf2[i][j] || fprintf(stderr, "%d %d %f %f\n", i, j, subm[i][j], mf2[i][j])<0);
		}
	}


}
