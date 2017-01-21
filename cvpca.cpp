#include <cxcore.h>
#include <cvext.h>

int main(int argc, char ** argv ){
	int m=5, n=10;
	cvRedirectError( cvSegReport );
	CvMat * data = cvCreateMat(m, n, CV_32F);
	CvMat * cov = cvCreateMat(n,n,CV_32F);
	CvMat * evect = cvCreateMat(MIN(m,n), n, CV_32F);
	CvMat * evect2 = cvCreateMat(m, m, CV_32F);
	CvMat * eval = cvCreateMat(MIN(m,n), 1, CV_32F);
	CvMat * mean = cvCreateMat(1, n, CV_32F);
	CvMat * AAT = cvCreateMat(m,m, CV_32F);
	CvMat * ATA = cvCreateMat(n,n, CV_32F);
	CvRNG rng = cvRNG();

	cvRandArr( &rng, data, CV_RAND_UNI, cvScalar(0), cvScalar(255) );
	cvCalcPCA( data, mean, eval, evect, CV_PCA_DATA_AS_ROW );
	cvCalcCovarMatrix( (const CvArr **) &data, 0, AAT, mean, CV_COVAR_ROWS | CV_COVAR_SCRAMBLED | CV_COVAR_SCALE);
	cvCalcCovarMatrix( (const CvArr **) &data, 0, ATA, mean, CV_COVAR_ROWS | CV_COVAR_SCALE);
	cvSaveMat5( "AAT.mat", AAT);
	cvSaveMat5( "ATA.mat", ATA);
	//cvSVD( AAT, eval, evect2, 0, CV_SVD_MODIFY_A | CV_SVD_U_T);
	cvSaveMat5( "A.mat", data );
	cvSaveMat5( "evect.mat", evect );
	cvSaveMat5( "mean.mat", mean );
	cvSaveMat5( "evals.mat", eval);
}
