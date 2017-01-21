#include <cxcore.h>
#include <cxmisc.h>
#include <cvtypes.h>
#include <matio.h>

CV_IMPL int cvSaveMat5(const char * fname, const CvArr * arr){
	int ret = 0;
	CV_FUNCNAME( "cvSaveMat5" );
	
	__BEGIN__;

	CvMat * mat = (CvMat *) arr;
	CvMat matStub;
	if( !CV_IS_MAT( mat ) ){
		mat = cvGetMat(mat, &matStub);
	}
	
	int array_class=0;
	int array_type=MAT_T_UNKNOWN;
	
	// dimensions sub-element
	int cn = CV_MAT_CN( mat->type );
	int width = mat->cols;
	int height = mat->rows;
	int step = mat->step;
	int dims[3] = {height, width, cn};
    int start[3]={0,0,0}, stride[3]={height,1,1}, edge[3]={1,width,cn};
    int ndim = cn > 1 ? 3 : 2;
	mat_t * matf;


	// determine datatype
	switch( CV_MAT_DEPTH(mat->type) ){
	case CV_8S:
		array_class = MAT_C_INT8;
		array_type = MAT_T_INT8;
		break;
	case CV_8U:
		array_class = MAT_C_UINT8;
		array_type = MAT_T_UINT8;
		break;
	case CV_16S:
		array_class = MAT_C_INT16;
		array_type = MAT_T_INT16;
		break;
	case CV_16U:
		array_class = MAT_C_UINT16;
		array_type = MAT_T_UINT16;
		break;
	case CV_32S:
		array_class = MAT_C_INT32;
		array_type = MAT_T_INT32;
		break;
	case CV_32F:
		array_class = MAT_C_SINGLE;
		array_type = MAT_T_SINGLE;
		break;
	case CV_64F:
		array_class = MAT_C_DOUBLE;
		array_type = MAT_T_DOUBLE;
		break;
	default:
		CV_ERROR( CV_StsUnsupportedFormat, "Could not determine type of matrix, it may be corrupt." );
	}

	// open file handle for writing
	matf = Mat_Open(fname,MAT_ACC_RDWR);
    if ( matf ) {
        matvar_t * matvar = Mat_VarCreate("cvarr",array_class,array_type,ndim,dims,NULL,0);
        if( Mat_VarWriteInfo( matf, matvar) != 0){
			CV_ERROR( CV_StsError, "Mat_VarWriteInfo failed");
		}
		for(int i=0; i<height; i++){
			start[0]=i;
			if( Mat_VarWriteData(matf, matvar, (void *)( mat->data.ptr+step*i ), start, stride, edge) != 0){
				CV_ERROR( CV_StsError, "Mat_VarWriteData failed");
			}
		}
        Mat_VarFree(matvar);
        Mat_Close(matf);
    } else {
		CV_ERROR( CV_StsUnsupportedFormat, "Error opening mat file for writing" );
    }

	

	__END__;
	return ret;
}

#ifdef MAIN
#include <highgui.h>
int main(int argc, char ** argv){
	IplImage * im = cvLoadImage(argv[1], 1);
	cvSaveMat5("A.mat", im);
}
#endif
