#include <cvext.h>
#include <cxcore.h>
#include <cv.h>
#include <cv.hpp>   // CvSepFilter
#include <stdio.h>
#include <stdarg.h>
#include <cxmisc.h>

#include <cv/cvext.hpp>
#include <cv/image.h>
#include <cv/matrix.hh>
#include <cv/cvarr_adapter.hpp>

#ifndef FILESYSTEM_PREFIX_LEN
# define FILESYSTEM_PREFIX_LEN(Filename) 0
#endif

#ifndef ISSLASH
# define ISSLASH(C) ((C) == '/')
#endif

CV_IMPL CvHashTable * cvCreateHashTable( CvTypeInfo * key_type, CvTypeInfo * value_type){
	return g_hash_table_new( NULL, NULL );
}
CV_IMPL int cvHashTableSize(CvHashTable * h){
	return g_hash_table_size( h );
}
CV_IMPL void * cvHashTableLookup(CvHashTable * h, const void * key){
	return g_hash_table_lookup(h, key);
}
CV_IMPL void cvHashTableInsert(CvHashTable * h, void * key, void * value){
	g_hash_table_insert(h, key, value);
}
CV_IMPL void cvReleaseHashTable( CvHashTable ** h ){
	g_hash_table_destroy( *h );
	*h = NULL;
}


CV_IMPL void cvGradient(const CvArr * srcarr, CvArr * dstarr_x, CvArr * dstarr_y, int aperture_size)
{
	CvSepFilter filter;
    
    CV_FUNCNAME( "cvGradient" );
    
    __BEGIN__;
    
    int origin = 0;
    int src_type, dst_x_type, dst_y_type;
    CvMat srcstub, *src = (CvMat*)srcarr;
    CvMat dststub_x, *dst_x = (CvMat*)dstarr_x;
    CvMat dststub_y, *dst_y = (CvMat*)dstarr_y;
    
    if( !CV_IS_MAT(src) )
        CV_CALL( src = cvGetMat( src, &srcstub ));
    if( CV_IS_IMAGE_HDR( srcarr ))
        origin = ((IplImage*)srcarr)->origin;

    src_type = CV_MAT_TYPE( src->type );

    if(dst_x){
		if( !CV_IS_MAT(dst_x) ) CV_CALL( dst_x = cvGetMat( dst_x, &dststub_x ));
    	
		dst_x_type = CV_MAT_TYPE( dst_x->type );
    	
		if( !CV_ARE_SIZES_EQ( src, dst_x ))
        	CV_ERROR( CV_StsBadArg, "src and dst_x have different sizes" );
		CV_CALL( filter.init_deriv( src->cols, src_type, dst_x_type, 1, 0,
				aperture_size, CvSepFilter::SUM_TO_1 | (origin ? CvSepFilter::FLIP_KERNEL : 0)));
		CV_CALL( filter.process( src, dst_x ));
	}

	if(dst_y){
	    if( !CV_IS_MAT(dst_y) ) CV_CALL( dst_y = cvGetMat( dst_y, &dststub_y ));
    	
		dst_y_type = CV_MAT_TYPE( dst_y->type );

		if( !CV_ARE_SIZES_EQ( src, dst_y ))
			CV_ERROR( CV_StsBadArg, "src and dst_y have different sizes" );

		CV_CALL( filter.init_deriv( src->cols, src_type, dst_y_type, 0, 1,
				aperture_size, CvSepFilter::SUM_TO_1 | (origin ? CvSepFilter::FLIP_KERNEL : 0)));
		CV_CALL( filter.process( src, dst_y ));
	}



	__END__;
}

CV_IMPL void cvImageSC(const char * win, CvArr * arr, double min, double max){
	IplImage * im = cvCreateImage(cvGetSize(arr), 8, 1);
	
	if(min==max){
		cvMinMaxLoc(arr, &min, &max);
	}

	if(min==max){
		cvZero(im);
	}
	else{
		cvScale(arr, im, 255.0/(max-min), 255*(-min)/(max-min));
	}
	cvShowImage(win, im);
	cvReleaseImage(&im);
}

CV_IMPL int cvIsEqual(CvArr * src_arr, CvArr * dst_arr){
	CvMat *src, *dst;
	CvMat src_stub;
	CvMat dst_stub;
	int rows, cols, step;
	uchar * dst_ptr;
	uchar * src_ptr;

	src = cvGetMat( src_arr, &src_stub );
	dst = cvGetMat( dst_arr, &dst_stub );

	if(src->rows!=dst->rows || src->cols!=dst->cols || src->type!=dst->type) return 0;

	rows = dst->rows;
	cols = dst->cols*CV_ELEM_SIZE(dst->type);
	step = dst->step;
	for(int i=0; i<rows; i++){
		dst_ptr = dst->data.ptr + step*i;
		src_ptr = src->data.ptr + step*i;
		for(int j=0; j<cols; j++){
			if( dst_ptr[j] != src_ptr[j] ) return 0;
		}
	}
	return 1;
}

// return number of array elements
// includes count of number of array channels
CV_IMPL int cvGetNumElements( CvArr * arr ){
	int size = -1;

	CV_FUNCNAME( "cvGetNumElements" );

	__BEGIN__;

	int type;
	int ndim;
	int sizes[CV_CN_MAX];

	CV_CALL( type = cvGetElemType( arr ) );
	CV_CALL( ndim = cvGetDims( arr, sizes ) );
	size=CV_MAT_CN(type);
	for(int i=0; i<ndim; i++){
		size*=sizes[i];
	}

	__END__;
	return size;
}

/** Examine header and if it matches rows,cols,type, return it,
 *  Otherwise, release header and call cvCreateMat(rows,cols,type) */
CV_IMPL CvMat * cvReCreateMat( CvMat * header, int rows, int cols, int type ){
	if(header && header->rows==rows && 
	   header->cols==cols && CV_MAT_TYPE(header->type)==type){
		return header;
	}
	if(!header){
		header = cvCreateMatHeader( rows, cols, type );
	}
	else{
		cvReleaseData( header );
		cvInitMatHeader( header, rows, cols, type );
	}
	cvCreateData( header );
	return header;
}

CV_IMPL void cvProject( CvArr * vec, CvArr * basis, CvArr * proj ){
	CvSize sz = cvGetSize(basis);
	CvMat *tmp = cvCreateMat( sz.height, 1, CV_32F );
	CvMat row;
	CvMat col;

	cvReshape(vec, &col, 1, sz.width);

	cvMatMul(basis, &col, tmp);

	cvReshape(tmp, &row, 1, 1);
	cvReshape(proj, &col, 1, 1);

	// 1xM MxN = 1xN
	cvMatMul(&row, basis, &col);
	
	cvReleaseMat( &tmp );
}

CV_IMPL void cvGramSchmidt( CvArr * basis, CvArr * orthoBasis ){
	CvSize sz = cvGetSize( basis );
	CvMat v, w;
	CvMat inc_basis;
	CvMat * tmp_vec = cvCreateMat( 1, sz.width, CV_32F );

	// w_1 = v_1 / ||v_1||
	cvGetRow( basis, &v, 0 );
	cvGetRow( orthoBasis, &w, 0);
	cvScale( &v, &w, 1/cvNorm(&v) );

	for(int i=1; i<sz.height; i++){
		// w_i* = v_i - proj(V_i-1, v_i)
		// w_i = w_i/||w_i||
		cvGetRow( basis, &v, i);
		cvGetRow( orthoBasis, &w, i);
		cvGetRows( orthoBasis, &inc_basis, 0, i );

		cvProject( &v, &inc_basis, tmp_vec );
		cvSub( &v, tmp_vec, &w );
		cvScale( &w, &w, 1/cvNorm(&w) );
	}
	
	cvReleaseMat( &tmp_vec );
}

#if 0

// for data -- each column is data -- 
// this is opposite the norm b/c of the way cvSVD calculates eigenvalues
// U can be MxN, but V must be NxN
// if computing EigenBasis of images NxN is intractable
//
// Computes principal component vectors v_i of matrix X
// s.t. X*v_i = c_i*v_i
//
//
CV_IMPL void cvCalcPCA( CvArr * data, CvArr * avg, CvArr * eigenvals, CvArr * eigenvects, int flags ){
	CV_FUNCNAME( "cvPCABasis" );

	CvMat col, col0, dataStub;
	CvMat * data0, *data_mat = (CvMat *) data;	
	CvMat avg_stub, *avg_mat = (CvMat *) avg;
	CvMat * S = NULL;
	bool data_as_col=(flags & CV_PCA_DATA_AS_COL);
	int numitems;
	CvScalar inv_numitems;

	__BEGIN__;

	if(! CV_IS_MAT( data ) ){
		data_mat = cvGetMat(data, &dataStub);
	}
	if(! CV_IS_MAT( avg_mat ) ){
		avg_mat = cvGetMat(avg, &avg_stub);
	}

	if(CV_MAT_DEPTH(avg_mat->type) != CV_32F && CV_MAT_DEPTH(avg_mat->type) != CV_64F){
		CV_ERROR(CV_BadDepth, "Mean must be CV_32F or CV_64F");
	}
	if(CV_MAT_DEPTH(data_mat->type) != CV_32F && CV_MAT_DEPTH(data_mat->type) != CV_64F){
		CV_ERROR(CV_BadDepth, "Data must be CV_32F or CV_64F");
	}

	// create array for scaled/shifted data 
	if(data_as_col){
		printf("%dx%dx%d\n", data_mat->rows, data_mat->cols, CV_MAT_CN(data_mat->type));
		if(CV_MAT_CN( data_mat->type ) > 1 ){
			CV_ERROR(CV_BadNumChannels, "Data in column format must be 1-channel");
		}
		data0 = cvCloneMat( (const CvMat *) data );
		printf("%dx%dx%d\n", data_mat->rows, data_mat->cols, CV_MAT_CN(data_mat->type));
	}
	else{
		if( CV_MAT_CN( data_mat->type ) > 1 ){
			data_mat = cvReshape( data, &dataStub, 1);
		}
		data0 = cvCreateMat( data_mat->cols, data_mat->rows, data_mat->type );
	}

	// fix avg to be 1-channel column vector
	if( CV_MAT_CN( avg_mat->type ) > 1 || avg_mat->cols > 1){
		avg_mat = cvReshape( avg, &avg_stub, 1, avg_mat->rows*avg_mat->cols*CV_MAT_CN(avg_mat->type));
	}

	if( avg_mat->rows != data0->rows ){
		printf("%dx%d %dx%d\n", data0->rows, data0->cols, avg_mat->rows, avg_mat->cols);
		CV_ERROR(CV_BadImageSize, "data dimension does not match mean dimensions");
	}

	if(!eigenvals){
		S = cvCreateMat( data0->cols, 1, data0->type );
		eigenvals = S;
	}

	
	// compute mean
	numitems = data0->cols;
	inv_numitems = cvScalarAll( 1.0/numitems );
	if((flags & CV_COVAR_USE_AVG)==0){
		cvZero(avg_mat);
		if(data_as_col){
			for(int i=0; i<numitems; i++){
				cvGetCol(data, &col, i);
				cvScaleAdd(&col, inv_numitems, avg_mat, avg_mat);
			}
		}
		else{
			CvMat row;
			for(int i=0; i<numitems; i++){
				cvGetRow(data, &row, i);
				cvReshape(&row, &col, 1, row.cols);
				cvScaleAdd(&col, inv_numitems, avg_mat, avg_mat);
			}
		}
	}
	cvCheckArr( avg_mat );

	// center data
	if(data_as_col){
		for(int i=0; i<numitems; i++){
			cvGetCol(data, &col, i);
			cvGetCol(data0, &col0, i);
			cvSub(&col, avg_mat, &col0);
		}
	}
	else{
		CvMat row;
		for(int i=0; i<numitems; i++){
			cvGetRow(data, &row, i);
			cvReshape(&row, &col, 1, data0->rows); // as a column vector
			cvGetCol(data0, &col0, i);
			cvSub(&col, avg_mat, &col0);
		}
	}

	cvCheckArr(data0);

	// compute eigen vectors
	//printf("U=%dx%d data=%dx%d\n", data0->rows, data0->cols, ((CvMat *)eigenvects)->rows, ((CvMat *)eigenvects)->cols);
	// compute A' = U * S * V'
	cvSVD( data0, eigenvals, eigenvects, NULL, CV_SVD_U_T + CV_SVD_MODIFY_A ); 

	cvReleaseMat( &data0 );

	if(S!=NULL){
		cvReleaseMat( &S ); 
	}
	else{
		if( (flags & CV_PCA_RET_SDV) != 0){
			// for standard deviation, need to scale
			cvScale(eigenvals, eigenvals, sqrt(1.0/(numitems-1)));
		}
		else{
			// singular values of A are square root of eigenvalues of A'A
			// for eigenvalues of covar matrix, need to scale by 1/nitems
			cvMul( eigenvals, eigenvals, eigenvals, 1.0/(numitems-1));
		}
	}

	__END__;
}
// return projection of data into eigenspace
CV_IMPL void cvProjectPCA( CvArr * data, CvArr * mean, CvArr * eigenVectors, CvArr * outResult, int dim)
{
	CV_FUNCNAME( "cvProjectPCA" );

	__BEGIN__;

    // y = A( x - mean )
    // where y = projection
    //       A = matrix of eigenVectors
    
    // check arguments
    CvMat dataMat;
    CvMat meanMat;
	CvMat truncEV;
	CvMat resMat;
    
    // reshape arrays into 1-channel row vectors 
	int data_size, mean_size, result_size;
	CV_CALL( data_size = cvGetNumElements( data ));
	CV_CALL( mean_size = cvGetNumElements( mean ));
	CV_CALL( result_size = cvGetNumElements( outResult ));

	if(mean_size!=data_size){
		CV_ERROR( CV_StsUnmatchedSizes, "data and mean do not have the same number of elements");
	}

    data = cvReshape( data, &dataMat, 1, data_size );
    mean = cvReshape( mean, &meanMat, 1, data_size );
    outResult = cvReshape( outResult, &resMat, 1, result_size);
    
	CvSize evSize;
	CV_CALL( evSize = cvGetSize( eigenVectors ) );
    // extract only top 'dim' eigenvectors
    if(dim>0 && dim<evSize.height){
        eigenVectors = cvGetSubRect( eigenVectors, &truncEV, cvRect(0,0,evSize.width, dim));
    }
        
    // Ax
    cvGEMM( eigenVectors, data, 1, NULL, 1, outResult);

    // -A*mean + Ax
    cvGEMM( eigenVectors, mean, -1, outResult, 1, outResult);

	__END__;
}
    
// return backprojection of point in eigenspace to reconstructed space
CV_IMPL
void cvBackProjectPCA( CvArr * eigenData, CvArr * mean, CvArr * eigenVectors, CvArr * outResult )
{
	CV_FUNCNAME( "cvEigenBackProject" );

	__BEGIN__;

	CvMat dataMat;
    CvMat meanMat;
    CvMat truncEV;
    CvMat resMat;
    

    // reshape matrices into 1-channel row vectors 
	int eigen_size, mean_size, result_size;
	CV_CALL( eigen_size = cvGetNumElements( eigenData ) );
	CV_CALL( mean_size = cvGetNumElements( mean ) );
	CV_CALL( result_size = cvGetNumElements( outResult ) );

	if(mean_size!=result_size){
		CV_ERROR( CV_StsUnmatchedSizes, "result and mean do not have the same number of elements");
	}

    eigenData = cvReshape( eigenData, &dataMat, 1, eigen_size ); 
    mean = cvReshape( mean, &meanMat, 1, result_size );
    outResult = cvReshape( outResult, &resMat, 1, result_size );
        
    // extract enough eigenvectors for the number of rows in eigenData
	CvSize edsz, evsz;
	CV_CALL( edsz = cvGetSize(eigenData) );
	CV_CALL( evsz = cvGetSize(eigenVectors) );
    if(edsz.height<evsz.height){
        eigenVectors = cvGetSubRect( eigenVectors, &truncEV, cvRect(0,0,evsz.width,edsz.height));
    }
    
    // x = A^T y + mean
    //std::cout<<"eigenVectors="<<*(Matrix<float> *)eigenVectors<<std::endl;
    //std::cout<<"eigenData="<<*(Matrix<float> *)eigenData<<std::endl;
    //std::cout<<"mean="<<*(Matrix<float> *)mean<<std::endl;
	//printf("cvEigenBackProject: coefs: %dx%d evect: %dx%d, mean: %dx%d result: %dx%d\n",
	//		dataMat.rows, dataMat.cols, truncEV.rows, truncEV.cols, meanMat.rows, meanMat.cols, resMat.rows, resMat.cols);
	
    //cvGEMM( eigenData, eigenVectors, 1, mean, 1, outResult, CV_GEMM_A_T | CV_GEMM_C_T );
    cvGEMM( eigenVectors, eigenData, 1, mean, 1, outResult, CV_GEMM_A_T );
    //std::cout<<"data="<<*(Matrix<float> *)outResult<<std::endl;
	
	__END__;
}
#endif

/** return only filename portion of path */
CVAPI(char *) cvBaseName(const char * name);
/** return only filename portion of path */
char * cvBaseName(const char * name){
	char const *base = name += FILESYSTEM_PREFIX_LEN (name);
	int all_slashes = 1;
	char const *p;

	for (p = name; *p; p++)
	{
		if (ISSLASH (*p))
			base = p + 1;
		else
			all_slashes = 0;
	}

	/* If NAME is all slashes, arrange to return `/'.  */
	if (*base == '\0' && ISSLASH (*name) && all_slashes)
		--base;

	return (char *) base;
}

/** extract filename without extension */
char * cvRootName(const char * name, char * dest){
	int len = strlen(name);
	strcpy( dest, name );
	for(int i=len-1; i>=0; i--){
		dest[i] = 0;
		if(name[i]=='.') break;
	}
	return dest;
}

/** return only directory portion of path */
char * cvDirName(const char * name, char * dest){
	char const *p;
	char const *end = name;
	char * ret = dest;
	
	// find last slash
	for(p=name; *p; p++){
		if(ISSLASH(*p)){
			end = p+1;
		}
	}
	// copy all characters up to there
	for(p=name; p!=end; p++){
		*dest=*p;
		dest++;
	}
	*dest='\0';
	return ret;
}

int icvInitExt(){
	cvRedirectError( cvSegReport );
	return 1;
}
int g_ext_is_init = icvInitExt();

int cvSegReport( int status, const char * func_name,
		const char *err_msg, const char * file_name,
		int line, void *userdata ){
	fprintf(stderr, "ERROR in %s(%s:%d) %s\n%s\n", func_name, file_name, line, cvErrorStr(status), err_msg);
	assert(0);
}

CvBox2D cvBox2DFromMoments( CvMoments * moments ){
	CvBox2D box;
	double m10,m00,m01;
	m00 = cvGetSpatialMoment(moments, 0, 0);
	m10 = cvGetSpatialMoment(moments, 1, 0);
	m01 = cvGetSpatialMoment(moments, 0, 1);

	double mu11 = moments->mu11;
	double mu20 = moments->mu20;
	double mu02 = moments->mu02;
	double inv_m00 = 1. / m00;
	box.center.x = m10 * inv_m00;
	box.center.y = m01 * inv_m00;
	double a = mu20 * inv_m00;
	double b = mu11 * inv_m00;
	double c = mu02 * inv_m00;

	/* Calculating width & height */
	double square = sqrt( 4 * b * b + (a - c) * (a - c) );

	/* Calculating orientation */
	box.angle = atan2( 2 * b, a - c + square );

	/* Calculating width & length of figure */
	double cs = cos( box.angle );
	double sn = sin( box.angle );

	double rotate_a = cs * cs * mu20 + 2 * cs * sn * mu11 + sn * sn * mu02;
	double rotate_c = sn * sn * mu20 - 2 * cs * sn * mu11 + cs * cs * mu02;
	box.size.width = sqrt( rotate_a * inv_m00 ) * 4;
	box.size.height = sqrt( rotate_c * inv_m00 ) * 4;
	return box;
}

void cvLocalMinMax(CvArr * arr, CvPoint p, CvSize win, double *min, double *max, CvPoint *minpos, CvPoint *maxpos){
	CvSize sz = cvGetSize( arr );
	int x0 = MAX(0, p.x-win.width/2);
	int x1 = MIN(p.x+win.width/2, sz.width);
	int y0 = MAX(0, p.y-win.height/2);
	int y1 = MIN(p.y+win.height/2, sz.height);
	CvRect roi = cvRect(x0, y0, x1-x0+1, y1-y0+1);
	CvMat mat;
	cvGetSubRect(arr, &mat, roi);
	cvMinMaxLoc( &mat, min, max, minpos, maxpos );
	if(minpos){
		minpos->x += roi.x;
		minpos->y += roi.y;
	}
	if(maxpos){
		maxpos->x += roi.x;
		maxpos->y += roi.y;
	}
}

IplImage ** cvCreateImagePyr( CvSize size, int depth, int nchannels, int steps ){
	IplImage ** pyr = (IplImage **) cvAlloc( sizeof(IplImage *)*(steps+1) );
	for(int i=0; i<steps; i++){
		pyr[i] = cvCreateImage( size, depth, nchannels );
		size.width /= 2;
		size.height /= 2;
		if(size.width==0 || size.height==0) break;
	}
	pyr[steps] = NULL;
	return pyr;
}

void cvReleaseImagePyr( IplImage *** _pyr ){
	IplImage ** pyr = *_pyr;
	for(int i=0; pyr[i]; i++){
		cvReleaseImage( &(pyr[i]) );
	}
	cvFree( _pyr );
}

void cvPyrUpAll( CvArr * src, CvArr ** dst, int filter ){
	// count number of levels
	int nlevels = 0;
	while(dst[0]!=NULL){
		nlevels++;
		dst++;
	}
	for(int i=nlevels-1; i>=0; i--){
		cvPyrUp(src, dst[i], filter);
		src=dst[i];
	}
}

void cvPyrDownAll( CvArr * src, CvArr ** dst, int filter ){
	CvSize src_sz, dst_sz;
	src_sz = cvGetSize(src);
	dst_sz = cvGetSize(dst[0]);
	if(src_sz.width == dst_sz.width){
		if(src!=dst[0]){
			cvCopy( src, dst[0] );
		}
		src=dst[0];
		dst++;
	}
	while(dst[0]!=NULL){
		cvPyrDown( src, dst[0], filter );
		src = dst[0];
		dst++;
	}
}

// Rearrange the quadrants of Fourier image so that the origin is at
// the image center
// src & dst arrays of equal size & type
void cvShiftDFT(CvArr * src_arr, CvArr * dst_arr ){
	CV_FUNCNAME("cvShiftDFT");
	
	__BEGIN__;

	CvMat q1stub, q2stub;
	CvMat q3stub, q4stub;
	CvMat * q1, * q2, * q3, * q4;

	CvSize size = cvGetSize(src_arr);
	CvSize dst_size = cvGetSize(dst_arr);
	int depth, cx, cy;

	if(dst_size.width != size.width || 
	   dst_size.height != size.height){
		CV_ERROR( CV_StsUnmatchedSizes, "" );	
	}

	cx = size.width/2;
	cy = size.height/2; // image center

	q1 = cvGetSubRect( src_arr, &q1stub, cvRect(0,0,size.width/2,size.height/2) );
	q2 = cvGetSubRect( src_arr, &q2stub, cvRect(cx,0,size.width/2,size.height/2) );
	q3 = cvGetSubRect( src_arr, &q3stub, cvRect(cx,cy,size.width/2,size.height/2) );
	q4 = cvGetSubRect( src_arr, &q4stub, cvRect(0,cy,size.width/2,size.height/2) );

	depth = CV_MAT_DEPTH(q1->type);


	if(src_arr==dst_arr){
		if( depth < CV_32F ){
			CV_ERROR( CV_StsUnsupportedFormat,
					"Only 32fC1, 64fC1 formats are supported" );
		}
		if( depth == CV_32F ){
			float tmp;
			SWAP_MAT_ELEMENTS( q1, q3, tmp );
			SWAP_MAT_ELEMENTS( q2, q4, tmp );
		}
		if( depth == CV_64F ){
			double tmp;
			SWAP_MAT_ELEMENTS( q1, q3, tmp );
			SWAP_MAT_ELEMENTS( q2, q4, tmp );
		}
	}
	else{
		CvMat dst;
		cvGetSubRect( dst_arr, &dst, cvRect(0,0,size.width/2,size.height/2) );
		if( !CV_ARE_TYPES_EQ( q3, &dst )){
			CV_ERROR( CV_StsUnmatchedFormats, "" );
		}
		cvCopy(q3, &dst);
		cvGetSubRect( dst_arr, &dst, cvRect(cx,0,size.width/2,size.height/2) );
		cvCopy(q4, &dst);
		cvGetSubRect( dst_arr, &dst, cvRect(cx,cy,size.width/2,size.height/2) );
		cvCopy(q1, &dst);
		cvGetSubRect( dst_arr, &dst, cvRect(0,cy,size.width/2,size.height/2) );
		cvCopy(q2, &dst);
	}

	__END__;
}

#if 0
/* Convert matrix to vector */
#define CV_MAT2VEC( mat, vdata, vstep, num )       \
    assert( (mat).rows == 1 || (mat).cols == 1 );  \
    (vdata) = ((mat).data.ptr);                    \
    if( (mat).rows == 1 )                          \
    {                                              \
        (vstep) = CV_ELEM_SIZE( (mat).type );      \
        (num) = (mat).cols;                        \
    }                                              \
    else                                           \
    {                                              \
        (vstep) = (mat).step;                      \
        (num) = (mat).rows;                        \
    }


#define ICV_RAND_SHUFFLE( suffix, type )                                                 \
static void icvRandShuffle_##suffix( uchar* data, size_t step, int num, CvRNG * rng )           \
{                                                                                        \
    type tmp;                                                                            \
    int i;                                                                               \
    int rn;                                                                              \
                                                                                         \
    for( i = 0; i < (num-1); i++ )                                                       \
    {                                                                                    \
        rn = cvRandInt( rng ) % (num - i);                                               \
        CV_SWAP( *((type*)(data + i * step)),                                            \
                 *((type*)(data + ( i + rn ) * step)),                                   \
                 tmp );                                                                  \
    }                                                                                    \
}   
        
ICV_RAND_SHUFFLE( 8U, uchar )
    
ICV_RAND_SHUFFLE( 16S, short )
        
ICV_RAND_SHUFFLE( 32S, int )

ICV_RAND_SHUFFLE( 32F, float )
    
void cvRandShuffle( CvMat* mat, CvRNG * rng)
{       
    CV_FUNCNAME( "cvRandShuffle" );

    __BEGIN__;

    uchar* data;
    size_t step;
    int num;
        
    if( (mat == NULL) || !CV_IS_MAT( mat ) || MIN( mat->rows, mat->cols ) != 1 )
    {       
        CV_ERROR( CV_StsUnsupportedFormat, "" );
    }       
            
    CV_MAT2VEC( *mat, data, step, num );
    switch( CV_MAT_TYPE( mat->type ) )
    {       
        case CV_8UC1:
            icvRandShuffle_8U( data, step, num, rng);
            break;
        case CV_16SC1:
            icvRandShuffle_16S( data, step, num, rng);
            break;
        case CV_32SC1:
            icvRandShuffle_32S( data, step, num, rng);
            break;
        case CV_32FC1:              
            icvRandShuffle_32F( data, step, num, rng);
            break;
        default: 
            CV_ERROR( CV_StsUnsupportedFormat, "" );
    }                           
    
    __END__;
} 
#endif

//	Returns the value ln[Gamma(xx)] for xx > 0.
double cvLogGamma(double xx)
{
//	    Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//			    accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double cvFactorial(int n)
{
    static int ntop=4;
    // Fill in table only as required.
    static float a[33]={1.0,1.0,2.0,6.0,24.0};
    int j;
    assert(n>=0);
    if (n > 32) return exp( cvLogGamma(n+1.0) );
    //Larger value than size of table is required. Actually, this big a value is going to overflow
    //on many computers, but no harm in trying.
	
	//Fill in table up to desired value.
    while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
    }
    return a[n];
}

//Returns the value ln( n! ) as a floating-point number.
double cvLogFactorial(int n)
{
  static double a[101];
  assert (n >= 0);
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]= cvLogGamma(n+1.0));
  else return cvLogGamma(n+1.0);
}

long cvBinomial(int n, int k){
	return cvFloor(0.5 + exp( cvLogFactorial(n) - cvLogFactorial(k) - cvLogFactorial(n-k)));
}

double cvCrossCorrelation(CvArr * src1, CvArr * src2, CvArr * tmp1, CvArr * tmp2){
	// xcorr = SUM[ (x(i) - m_x) * (y(i)-m_y) ] / sqrt( SUM[ (x(i)-m_x)^2 ] * SUM[ (y(i)-m_y)^2 ] )
	
	// m_x, m_y 
	CvScalar mean1 = cvAvg( src1 );
	CvScalar mean2 = cvAvg( src2 );
	
	// x(i) - m_x 
	cvSubS(src1, mean1, tmp1);
	cvSubS(src2, mean2, tmp2);
	
	//  sqrt( SUM[ (x(i)-m_x)^2 ] * SUM[ (y(i)-m_y)^2 ] )
	double bottom = cvNorm( tmp1, NULL, CV_L2 ) * cvNorm( tmp2, NULL, CV_L2 );

	//  (x(i) - m_x) * (y(i)-m_y)
	cvMul( tmp1, tmp2, tmp1 );

	// SUM[ (x(i) - m_x) * (y(i)-m_y) ]
	CvScalar top = cvSum( tmp1 );

	return top.val[0]/bottom;
}
struct cvIsLocalPeak2Adapter {
	bool is_min;  // flag to toggle min or max peak
	bool test_center;
	int size;
	IplImage * arr_cmp;
	IplImage * res;
	cv::Image<uchar, 1> * mask;
	cvIsLocalPeak2Adapter():
		test_center(false),
		size(3),
		arr_cmp(NULL),
		res(NULL),
		mask(NULL)
	{
	}
	template <typename T, int nch>
	int operator() (const cv::Image<T, nch> & im ){
		if(is_min){
			if(this->arr_cmp){
				cvIsLocalMin2( im, *cv::Image<T,nch>::safe_cast(this->arr_cmp), 
					       *cv::Image<uchar, 1>::safe_cast(this->res), 
						   this->test_center, this->size, this->mask );
			}
			else{
				cvIsLocalMin( im, *cv::Image<uchar, 1>::safe_cast(this->res), this->size, this->mask );
			}
		}
		else{
			if(this->arr_cmp){
				cvIsLocalMax2( im, *cv::Image<T,nch>::safe_cast(this->arr_cmp), 
						       *cv::Image<uchar, 1>::safe_cast(this->res),
							   this->test_center, this->size, this->mask );
			}
			else{
				cvIsLocalMax( im, *cv::Image<uchar, 1>::safe_cast(this->res), this->size, this->mask);
			}
		}
		return 0;
	}
};

CV_IMPL void cvIsLocalMin2( CvArr * arr, CvArr * arr_cmp, CvArr * dst, int test_center, int size, CvArr * mask) 
{

	IplImageAdapter1< cvIsLocalPeak2Adapter > adapter;
	IplImage arr_cmp_stub;
	IplImage dst_stub;
	IplImage mask_stub;
	
	adapter.func.arr_cmp =  cvGetImage( arr_cmp, &arr_cmp_stub );
	adapter.func.res = cvGetImage( dst, &dst_stub );
	adapter.func.test_center = (test_center!=0);
	adapter.func.size = size;
	adapter.func.mask = (mask ? cv::Image<uchar, 1>::safe_cast( cvGetImage( mask, &mask_stub ) ): NULL);
	adapter.func.is_min = true;

	adapter( arr );
}

CV_IMPL void cvIsLocalMax2( CvArr * arr, CvArr * arr_cmp, CvArr * dst, int test_center, int size, CvArr * mask) 
{
	IplImageAdapter1< cvIsLocalPeak2Adapter > adapter;
	IplImage arr_cmp_stub;
	IplImage dst_stub;
	IplImage mask_stub;
	
	adapter.func.arr_cmp =  cvGetImage( arr_cmp, &arr_cmp_stub );
	adapter.func.res = cvGetImage( dst, &dst_stub );
	adapter.func.test_center = (test_center!=0);
	adapter.func.size = size;
	adapter.func.mask = (mask ? cv::Image<uchar, 1>::safe_cast( cvGetImage( mask, &mask_stub ) ): NULL);
	adapter.func.is_min = false;

	adapter( arr );
}

CV_IMPL void cvIsLocalMin( CvArr * arr, CvArr * dst, int size, CvArr * mask )
{
	IplImageAdapter1< cvIsLocalPeak2Adapter > adapter;
	IplImage dst_stub;
	IplImage mask_stub;

	adapter.func.res = cvGetImage( dst, &dst_stub );
	adapter.func.size = size;
	adapter.func.mask = (mask ? cv::Image<uchar, 1>::safe_cast( cvGetImage( mask, &mask_stub ) ): NULL);
	adapter.func.is_min = true;

	adapter( arr );
}

CV_IMPL void cvIsLocalMax( CvArr * arr, CvArr * dst, int size, CvArr * mask )
{
	IplImageAdapter1< cvIsLocalPeak2Adapter > adapter;
	IplImage dst_stub;
	IplImage mask_stub;

	adapter.func.res = cvGetImage( dst, &dst_stub );
	adapter.func.size = size;
	adapter.func.mask = (mask ? cv::Image<uchar, 1>::safe_cast( cvGetImage( mask, &mask_stub ) ): NULL);
	adapter.func.is_min = false;

	adapter( arr );
}

struct cvMinMaxTopNAdapter {
	int num;
	float min; 
	float max;
	CvPoint min_point;
	CvPoint max_point;
	bool invert;

	cvMinMaxTopNAdapter() :
		num(50),
		invert(false)
	{}

	template<typename T, int nch>
	int operator() ( const cv::Image<T, nch> & im ){
		cvMinMaxTopN( im, num, min, max, min_point, max_point, invert );
		return 0;
	}
};

void cvMinMaxTopN( CvArr * arr, int num, float * min, float * max, CvPoint * minp, CvPoint * maxp, int invert){
	IplImageAdapter1< cvMinMaxTopNAdapter > func;
	func.func.num = num;
	func.func.invert = (invert != 0);
	func( arr );
	if(min) *min = func.func.min;
	if(max) *max = func.func.max;
	if(minp) *minp = func.func.min_point;
	if(maxp) *maxp = func.func.max_point;
}

// solve for matrix L s.t. A = L*L^T  where L is lower triangular and L^T is upper triangular
// formula is as follows
//
// L(i,j) = 1/L(i,j) * ( A(i,j) - SUM(k=0...j-1, L(i,k)*L(j,k)) )  for i > j
// L(i,i) = sqrt( A(i,i) - SUM(k=0...i-1, L(i,k)*L(i,k)) ) for i==j

// CvArr implementation
CVARR_TO_MAT_FUNC_BINARY( cvCholesky );

void cvRotateZoomImage(const CvArr * src, CvArr * dest, CvPoint2D32f center, float angle, float scale){
	// rotation matrix
	Matrix<float> r(3,3);
	Matrix<float> ri(3,3);
	Matrix<float> M(2,3);

	r.identity();
	r.rotateObject(2,angle+M_PI);
	ri = r.inv();
	scale = 1.0/scale;

	M(0,0) = -ri(0,0)*scale; M(0,1) = ri(0,1)*scale; M(0,2) = center.x;
	M(1,0) = ri(1,0)*scale; M(1,1) = -ri(1,1)*scale; M(1,2) = center.y;

	cvGetQuadrangleSubPix(src, dest, &M);
}
void cvRotateZoomImage2(const CvArr * src, CvArr * dest, CvPoint2D32f center, float angle, float scale, int flags, CvScalar fillval){
	CvSize dsize = cvGetSize(dest);
	Matrix<float> M(2,3);
	float a,b;

	a = scale*cos(angle);
	b = scale*sin(angle);
	M(0,0) =  a; M(0,1) = b; M(0,2) = dsize.width/2 - b*center.y - a*center.x;
	M(1,0) = -b; M(1,1) = a; M(1,2) = dsize.height/2 + b*center.x - a*center.y;
	
	cvWarpAffine( src, dest, &M, flags, fillval );
}

#if 0
CvMat * cvWarpAffineQMatrix( const CvPoint2D32f * src, const CvPoint2D32f * dst, CvMat * map_matrix ){
	// Want to solve for the matrix A s.t. 
    // for all points X=[x y] in src and points X'=[x', y'] in dst
    // A * ~X = X'
    // where ~X is the augmented vector [ x y 1 ]
    // and
    //  A = [ a b c ]
    //      [ d e f ]

    // for each point
    // A * ~X = X' gives two constraints,
    // ax + by + c = x'
    // dx + ey + f = y'

    // by rearranging terms and stacking constraints we get a linear system of equations of the form
    // H * ~A = B
    //
    // where 
    // ~A = [ a b c d e f ]^T
    //  B = [ x' y' ... ]^T  for each point
    //  H = [ x y 1 0 0 0 ]
    //      [ 0 0 0 x y 1 ]
    //      [ ..........  ]  for each point

    // with three or more constraints, ~A can be solved for directly.

    CV_FUNCNAME( "cvWarpAffineQMatrix" );

    __BEGIN__;

    CvMat mA, mX, mB;
    double A[6*6];
    double B[6];
    double x[6];
	
    cvInitMatHeader(&mA, 6, 6, CV_64F, A);
    cvInitMatHeader(&mB, 6, 1, CV_64F, B);
    cvInitMatHeader(&mX, 6, 1, CV_64F, x);

	if( !src || !dst || !map_matrix )
		CV_ERROR( CV_StsNullPtr, "" );

    for(int i=0; i<3; i++){
        int j = i*12;
        int k = i*12+6;
        A[j] = A[k+3] = src[i].x;
        A[j+1] = A[k+4] = src[i].y;
        A[j+2] = A[k+5] = 1;
        A[j+3] = A[j+4] = A[j+5] = 0;
        A[k] = A[k+1] = A[k+2] = 0;
        B[i*2] = dst[i].x;
        B[i*2+1] = dst[i].y;
    }
    cvSolve(&mA, &mB, &mX);

    mX = cvMat( 2, 3, CV_64FC1, x );
    cvConvert( &mX, map_matrix );

    __END__;
    return map_matrix;
}
#endif

/*Function///////////////////////////////////////////////////////////////

Name:       cvShowManyImages

Purpose:    This is a function illustrating how to display more than one 
               image in a single window using Intel OpenCV

Parameters: char *title: Title of the window to be displayed
            int nArgs:   Number of images to be displayed
            ...:         IplImage*, which contains the images

Language:   C++

The method used is to set the ROIs of a Single Big image and then resizing 
and copying the input images on to the Single Big Image.

This function does not stretch the image... 
It resizes the image without modifying the width/height ratio..

This function can be called like this:

cvShowManyImages("Images", 2, img1, img2);
or
cvShowManyImages("Images", 5, img2, img2, img3, img4, img5);

This function can display upto 12 images in a single window.
It does not check whether the arguments are of type IplImage* or not.
The maximum window size is 700 by 660 pixels.
Does not display anything if the number of arguments is less than
    one or greater than 12.

If you pass a pointer that is not IplImage*, Error will occur.
Take care of the number of arguments you pass, and the type of arguments, 
which should be of type IplImage* ONLY.

Idea was from BettySanchi of OpenCV Yahoo! Groups.

If you have trouble compiling and/or executing
this code, I would like to hear about it.

You could try posting on the OpenCV Yahoo! Groups
[url]http://groups.yahoo.com/group/OpenCV/messages/ [/url]


Parameswaran, 
Chennai, India.

cegparamesh[at]gmail[dot]com            

...
///////////////////////////////////////////////////////////////////////*/

CV_IMPL void cvShowManyImages(char* title, CvArr ** images, int nArgs) {

    // img - Used for getting the arguments 
    IplImage *img;
	IplImage header;

    // DispImage - the image in which input images are to be copied
    IplImage *DispImage;

    int size;
    int i;
    int m, n;
    int x, y;

    // w - Maximum number of images in a row 
    // h - Maximum number of images in a column 
    int w, h;

    // scale - How much we have to resize the image
    float scale;
    int max;

    // If the number of arguments is lesser than 0 or greater than 12
    // return without displaying 
    if(nArgs <= 0) {
        printf("Number of arguments too small....\n");
        return;
    }
    else if(nArgs > 12) {
        printf("Number of arguments too large....\n");
        return;
    }
    // Determine the size of the image, 
    // and the number of rows/cols 
    // from number of arguments 
    else if (nArgs == 1) {
        w = h = 1;
        size = 300;
    }
    else if (nArgs == 2) {
        w = 2; h = 1;
        size = 300;
    }
    else if (nArgs == 3 || nArgs == 4) {
        w = 2; h = 2;
        size = 300;
    }
    else if (nArgs == 5 || nArgs == 6) {
        w = 3; h = 2;
        size = 200;
    }
    else if (nArgs == 7 || nArgs == 8) {
        w = 4; h = 2;
        size = 200;
    }
    else {
        w = 4; h = 3;
        size = 150;
    }

    // Create a new 3 channel image
    DispImage = cvCreateImage( cvSize(100 + size*w, 60 + size*h), 8, 3 );

    // Loop for nArgs number of arguments
    for (i = 0, m = 20, n = 20; i < nArgs; i++, m += (20 + size)) {

        // Get the Pointer to the IplImage
        img = (IplImage *)images[i]; 
		
		if(!CV_IS_IMAGE(img)){
			img = cvGetImage(images[i], &header);
		}

        // Check whether it is NULL or not
        // If it is NULL, release the image, and return
        if(img == 0) {
            printf("Invalid arguments");
            cvReleaseImage(&DispImage);
            return;
        }

		// convert image if not 8UC3
		if(cvGetElemType(img)!=CV_8UC3){
			img = cvCreateImage( cvGetSize(img), 8, 3 );
			cvConvertImage( images[i], img );
		}

        // Find the width and height of the image
        x = img->width;
        y = img->height;

        // Find whether height or width is greater in order to resize the image
        max = (x > y)? x: y;

        // Find the scaling factor to resize the image
        scale = (float) ( (float) max / size );

        // Used to Align the images
        if( i % w == 0 && m!= 20) {
            m = 20;
            n+= 20 + size;
        }

        // Set the image ROI to display the current image
        cvSetImageROI(DispImage, cvRect(m, n, (int)( x/scale ), (int)( y/scale )));

        // Resize the input image and copy the it to the Single Big Image
        cvResize(img, DispImage);

        // Reset the ROI in order to display the next image
        cvResetImageROI(DispImage);

		if( img!=images[i] && img!=&header ){
			cvReleaseImage( &img );
		}
    }

    // Create a new window, and show the Single Big Image
    cvShowImage( title, DispImage);

    // Release the Image Memory
    cvReleaseImage(&DispImage);
}

#if 0
CV_IMPL CvArr * cvRange(CvArr * vec, double start, double end){
	CvMat mat_stub, *mat=(CvMat *)vec;
	double step;
	float *ptr;
	double val=start;
	int i,j;
	int rows,cols;
	
	if(!CV_IS_MAT(mat)) mat = cvGetMat(vec, &mat_stub);

	rows = mat->rows;
	cols = mat->cols;
	step = (end-start)/mat->cols;
	for(i=0; i<rows; i++){
		ptr = (float *)(mat->data.ptr + mat->step*i);
		for(j=0; j<cols; j++){
			*ptr = val;
			val+=step;
			++ptr;
		}
	}
	return vec;
}
#endif

template <typename T>
void icvSortInPlace( T * data, int len ){
	T tmp;
	/*
	std::cout<<"Original:";
	for(int ii=0; ii<len; ii++){
		std::cout<<" "<<data[ii];
	}
	std::cout<<std::endl;
	*/
	for(int i=0; i<len; i++){
		for(int j=i-1; j>=0; j--){
			if(data[j+1]<data[j]) CV_SWAP(data[j+1],data[j],tmp);
		}
	}
	/*std::cout<<"Sorted:";
	for(int ii=0; ii<len; ii++){
		std::cout<<" "<<data[ii];
	}
	std::cout<<std::endl;
	*/
}

#if 0
/**
 * Pulled from http://en.wikipedia.org/wiki/Shell_sort
 * BY BSF
 * Performs an insertion sort on elements of a[] with the given gap.
 * If gap == 1, performs an ordinary insertion sort.
 * If gap >= length, does nothing.
 * Based on a C implementation from the article Insertion sort implementations
 */
template <typename T>
void icvShellSortPhase(T a[], int length, int gap) {
    int i;
    for (i = gap; i < length; ++i) {
        T value = a[i];
        int j;
        for (j = i - gap; j >= 0 && a[j] > value; j -= gap) {
            a[j + gap] = a[j];
        }
        a[j + gap] = value;
    }
}
    
template <typename T>
void icvShellSort(T a[], size_t length) {
    /*
     * gaps[] should approximate a geometric progression.
     * The following sequence is the best known in terms of
     * the average number of key comparisons made [2]
     */
    static const int gaps[] = {
        1, 4, 10, 23, 57, 132, 301, 701, 0
    };
    int sizeIndex;
    
    for (sizeIndex = sizeof(gaps)/sizeof(gaps[0]) - 1;
               sizeIndex >= 0;
               --sizeIndex)
        icvShellSortPhase(a, length, gaps[sizeIndex]);
}

/// modification of basic shell sort to modify array of indices instead of actual array
/// step is in # of integers
template <typename T>
void icvShellSortIdx(const T a[], int * idx, int step, size_t length){
	static const int gaps[] = {
		701, 301, 132, 57, 23, 10, 4, 1, 0 
	};
	for(int sizeIndex=0; gaps[sizeIndex]!=0; sizeIndex++){
		for(int i=gaps[sizeIndex]; i < length; ++i){
			int vidx = idx[i*step];
			for(int j = i-gap; j>=0 && a[ idx[j] ] > a[ vidx ]; j-=gap){
				idx[ (j+gap)*step ] = idx[j*step];
			}
			idx[ (j + gap)*step ] = vidx;
		}
	}
}

void cvSortIdx(CvArr * arr, CvArr * idx, int flags){
	CvMat vec_stub, *vec;
	CvMat idx_stub, *idx_vec;

	int type = cvGetElemType(arr);
	CvSize sz = cvGetSize(arr);
	uchar * ptr;
	uchar * iptr;
	int istep;
	int step;
	int elem_size = CV_ELEM_SIZE( type );

	vec = cvReshape( arr, &vec_stub, 0, sz.width*sz.height );
	idx_vec = cvReshape( idx, &idx_stub, 0, sz.width*sz.height );

	cvGetRawData( vec, &ptr, &step );
	cvGetRawData( idx, &iptr, &istep );
	istep/=sizeof(int);
	step/=elem_size;

	if(flags==0){
		for(int y=0; y<size.height; y++){
			int * irow = (int *)( iptr + y*istep );
			for(int x=0; x<size.width; x++){
				irow[x] = y*step+x;
			}
		}
	}

	switch(type){
		case CV_8U:
			icvShellSortIdx( (uchar *)ptr, (int *)iptr, istep, sz.width*sz.height );
			break;
		case CV_8S:
			icvShellSortIdx( (char *)ptr, (int *)iptr, istep, sz.width*sz.height );
			break;
		case CV_16S:
			icvShellSortIdx( (short *)ptr, (int *)iptr, istep, sz.width*sz.height );
			break;
		case CV_32F:
			icvShellSortIdx( (float *)ptr, (int *)iptr, istep, sz.width*sz.height );
			break;
		case CV_32S:
			icvShellSortIdx( (int *)ptr, (int *)iptr, istep, sz.width*sz.height );
			break;
		default:
			break;
	}
}
#endif

#if 0
template <typename T>
void icvSort5( T * data, int nch, int len, T * sorted ){
	assert(len <= 5);
	assert(len > 0 );
	for(int i=0; i<len; i++){
		sorted[i] = data[i*nch];
		// check old elements
		for(int j=0; j<i; j++){
			if(data[i*nch] < sorted[j]){
				// shift all elements
				for(int k=i; k>j; k--){
					sorted[k]=sorted[k-1];
				}
				// insert element
				sorted[j] = data[i*nch];

				// terminate j loop
				j=i;
			}
		}
	}
	/*std::cout<<"Original:";
	for(int ii=0; ii<len; ii++){
		std::cout<<" "<<data[ii*nch];
	}
	std::cout<<std::endl;
	std::cout<<"Sorted:";
	for(int ii=0; ii<len; ii++){
		std::cout<<" "<<sorted[ii];
	}
	std::cout<<std::endl;*/
}
#endif

template <typename T>
T icvSelect( T * data, int nch, int rows, int cols, int step, int n ){
	int i,j;
	int extra = (rows*cols)%5 != 0 ? 1 : 0;
	T sorted[5];
	//memset(sorted, 0, sizeof(T)*5);


	// base case 
	if( rows*cols <= 5 ){
		int k=0;
		for(int i=0; i<rows; i++){
			T * ptr = (T*)((uchar *)data+i*step);
			for(int j=0; j<cols; j++){
				sorted[k++]=ptr[j*nch];
			}
		}
		icvSortInPlace( sorted, k );
		return sorted[ n ];
	}

	int nmedians = rows*cols/5+extra;
	T * median5 = (T*) cvAlloc( sizeof(T)*(nmedians) );

	// 1. split into len/5 sub arrays of length 5 and one of n mod 5
	//    and find the median of each
	int k=0;
	for(i=0; i<rows; i++){
		T * ptr = (T*)((uchar *)data+i*step);
		for(j=0; j<cols; j++){
			sorted[k%5] = ptr[j*nch];	
			++k;
			if(k%5==0){
				icvSortInPlace( sorted, 5 );
				median5[ (k/5)-1 ] = sorted[ 2 ];
			}
		}
	}
	if(k%5!=0){
		int nelem = k%5;
		icvSortInPlace( sorted, nelem );
		median5[ (k/5) ] = sorted[ nelem/2 ];
	}

	T median = icvSelect( median5, 1, nmedians, 1, sizeof(T), nmedians/2 );

	// partition input array around median of medians and recurse on
	// partition containing Nth element
	T * part1 = (T*) cvAlloc( sizeof(T)*rows*cols );
	T * part2 = (T*) cvAlloc( sizeof(T)*rows*cols );
	int i1=0,i2=0;
	for(i=0; i<rows; i++){
		T * ptr = (T*)((uchar *)data+i*step);
		for(j=0; j<cols; j++){
			if(ptr[j*nch] < median){
				part1[i1++]=ptr[j*nch];
			}
			else{
				part2[i2++]=ptr[j*nch];
			}
		}
	}

	T selected;
	if( i1 == 0 ){
		selected = median;
	}
	else if( i2 == 0 ){
		selected = median;
	}
	else if( n < i1 ){
		selected = icvSelect( part1, 1, i1, 1, sizeof(T), n );
	}
	else{
		selected = icvSelect( part2, 1, i2, 1, sizeof(T), n-i1 );
	}

	cvFree( &median5 );
	cvFree( &part1 );
	cvFree( &part2 );

	return selected;
}

CV_IMPL CvScalar cvSelect( const CvArr * arr, int n ){
	CvMat mat_stub, *mat = (CvMat *)arr;
	if(!CV_IS_MAT(mat)) mat = cvGetMat(arr, &mat_stub);

	CvScalar selected;
	int nch = CV_MAT_CN(mat->type);
	int k;

	CV_FUNCNAME("cvSelect");
	__BEGIN__;
	for( k=0; k<nch; k++){
		switch(CV_MAT_DEPTH(mat->type)){
			case CV_8U:
				selected.val[k] = icvSelect( mat->data.ptr+k, nch, mat->rows, mat->cols, mat->step, n );
				break;
			case CV_8S:
				selected.val[k] = icvSelect( (char *)(mat->data.ptr+k), nch, mat->rows, mat->cols, mat->step, n );
				break;
			case CV_16S:
				selected.val[k] = icvSelect( mat->data.s+k, nch, mat->rows, mat->cols, mat->step, n );
				break;
			case CV_32S:
				selected.val[k] = icvSelect( mat->data.i+k, nch, mat->rows, mat->cols, mat->step, n );
				break;
			case CV_32F:
				selected.val[k] = icvSelect( mat->data.fl+k, nch, mat->rows, mat->cols, mat->step, n );
				break;
			case CV_64F:
				selected.val[k] = icvSelect( mat->data.ptr+k, nch, mat->rows, mat->cols, mat->step, n );
				break;
			default:
				CV_ERROR( CV_BadDepth, "" );
		}
	}
	__END__;

	return selected;
}

CV_IMPL CvScalar cvMedian( const CvArr * arr ){
	CvSize sz = cvGetSize(arr);
	return cvSelect( arr, sz.width*sz.height/2 ); 
}


CV_IMPL void cvReducePtr( const CvArr * srcarr, CvArr * dstarr, CvReduceFunc func, int dim ){
	CvMat src_stub, *src=(CvMat *) srcarr;
	CvMat dst_stub, *dst=(CvMat *) dstarr;

	CV_FUNCNAME("cvReduce");
	__BEGIN__;

	if(!CV_IS_MAT(src)) src = cvGetMat( srcarr, &src_stub );

	if(!CV_IS_MAT(dst)) dst = cvGetMat( dstarr, &dst_stub );

	if(dim>=2){
		CV_ERROR( CV_StsUnsupportedFormat, "dim >= 2 not supported" );
	}

	dst = cvReshape( dstarr, &dst_stub, 0, 1 );

	if(dim==0){
		CvMat row;
		if( src->rows != dst->cols ){
			CV_ERROR( CV_StsBadSize, "number of elements of dstarr must equal # of rows of src\n" );
		}
		for(int i=0; i<src->rows; i++){
			cvGetRow( src, &row, i);
			cvSet1D( dst, i, func( &row ) );
		}
	}
	else if(dim==1){
		CvMat col;
		if( src->cols != dst->cols ){
			CV_ERROR( CV_StsBadSize, "number of elements of dstarr must equal # of rows of src\n" );
		}
		for(int i=0; i<src->cols; i++){
			cvGetCol( src, &col, i);
			cvSet1D( dst, i, func( &col ) );
		}
	}

	__END__;
}

CV_IMPL void
cvKMeans3( const CvArr* samples_arr, int cluster_count,
           CvArr * centers_arr, CvArr* labels_arr, CvTermCriteria termcrit )
{
    CvMat* centers = 0;
    CvMat* old_centers = 0;
    CvMat* counters = 0;

    CV_FUNCNAME( "cvKMeans2" );

    __BEGIN__;

    CvMat samples_stub, labels_stub;
    CvMat* samples = (CvMat*)samples_arr;
    CvMat* labels = (CvMat*)labels_arr;
    CvMat* temp = 0;
    CvRNG rng = CvRNG(-1);
    int i, j, k, sample_count, dims;
    int ids_delta, iter;
    double max_dist;

    if( !CV_IS_MAT( samples ))
        CV_CALL( samples = cvGetMat( samples, &samples_stub ));

    if( !CV_IS_MAT( labels ))
        CV_CALL( labels = cvGetMat( labels, &labels_stub ));

    if( cluster_count < 1 )
        CV_ERROR( CV_StsOutOfRange, "Number of clusters should be positive" );

    if( CV_MAT_DEPTH(samples->type) != CV_32F || CV_MAT_TYPE(labels->type) != CV_32SC1 )
        CV_ERROR( CV_StsUnsupportedFormat,
        "samples should be floating-point matrix, cluster_idx - integer vector" );

    if( labels->rows != 1 && (labels->cols != 1 || !CV_IS_MAT_CONT(labels->type)) ||
        labels->rows + labels->cols - 1 != samples->rows )
        CV_ERROR( CV_StsUnmatchedSizes,
        "cluster_idx should be 1D vector of the same number of elements as samples' number of rows" );

    CV_CALL( termcrit = cvCheckTermCriteria( termcrit, 1e-6, 100 ));

    termcrit.epsilon *= termcrit.epsilon;
    sample_count = samples->rows;

    if( cluster_count > sample_count )
        cluster_count = sample_count;

    dims = samples->cols*CV_MAT_CN(samples->type);
    ids_delta = labels->step ? labels->step/(int)sizeof(int) : 1;

    CV_CALL( centers = cvCreateMat( cluster_count, dims, CV_64FC1 ));
    CV_CALL( old_centers = cvCreateMat( cluster_count, dims, CV_64FC1 ));
    CV_CALL( counters = cvCreateMat( 1, cluster_count, CV_32SC1 ));

	cvConvert( centers_arr, centers );

    counters->cols = cluster_count; // cut down counters
    max_dist = termcrit.epsilon*2;

    for( iter = 0; iter < termcrit.max_iter; iter++ )
    {
        // assign labels
        for( i = 0; i < sample_count; i++ )
        {
            float* s = (float*)(samples->data.ptr + i*samples->step);
            int k_best = 0;
            double min_dist = DBL_MAX;

            for( k = 0; k < cluster_count; k++ )
            {
                double* c = (double*)(centers->data.ptr + k*centers->step);
                double dist = 0;
                
                j = 0;
                for( ; j <= dims - 4; j += 4 )
                {
                    double t0 = c[j] - s[j];
                    double t1 = c[j+1] - s[j+1];
                    dist += t0*t0 + t1*t1;
                    t0 = c[j+2] - s[j+2];
                    t1 = c[j+3] - s[j+3];
                    dist += t0*t0 + t1*t1;
                }

                for( ; j < dims; j++ )
                {
                    double t = c[j] - s[j];
                    dist += t*t;
                }

                if( min_dist > dist )
                {
                    min_dist = dist;
                    k_best = k;
                }
            }

            labels->data.i[i*ids_delta] = k_best;
        }

        if( max_dist < termcrit.epsilon )
            break;

        CV_SWAP( centers, old_centers, temp );
        
		// recompute centers
        cvZero( centers );
        cvZero( counters );

        for( i = 0; i < sample_count; i++ )
        {
            float* s = (float*)(samples->data.ptr + i*samples->step);
            k = labels->data.i[i*ids_delta];
            double* c = (double*)(centers->data.ptr + k*centers->step);
            for( j = 0; j <= dims - 4; j += 4 )
            {
                double t0 = c[j] + s[j];
                double t1 = c[j+1] + s[j+1];

                c[j] = t0;
                c[j+1] = t1;

                t0 = c[j+2] + s[j+2];
                t1 = c[j+3] + s[j+3];

                c[j+2] = t0;
                c[j+3] = t1;
            }
            for( ; j < dims; j++ )
                c[j] += s[j];
            counters->data.i[k]++;
        }

        if( iter > 0 )
            max_dist = 0;

        for( k = 0; k < cluster_count; k++ )
        {
            double* c = (double*)(centers->data.ptr + k*centers->step);
            if( counters->data.i[k] != 0 )
            {
                double scale = 1./counters->data.i[k];
                for( j = 0; j < dims; j++ )
                    c[j] *= scale;
            }
            else
            {
                i = cvRandInt( &rng ) % sample_count;
                float* s = (float*)(samples->data.ptr + i*samples->step);
                for( j = 0; j < dims; j++ )
                    c[j] = s[j];
            }
            
            if( iter > 0 )
            {
                double dist = 0;
                double* c_o = (double*)(old_centers->data.ptr + k*old_centers->step);
                for( j = 0; j < dims; j++ )
                {
                    double t = c[j] - c_o[j];
                    dist += t*t;
                }
                if( max_dist < dist )
                    max_dist = dist;
            }
        }


    }

    cvZero( counters );
    for( i = 0; i < sample_count; i++ )
        counters->data.i[labels->data.i[i]]++;

    // ensure that we do not have empty clusters
    for( k = 0; k < cluster_count; k++ )
        if( counters->data.i[k] == 0 )
            for(;;)
            {
                i = cvRandInt(&rng) % sample_count;
                j = labels->data.i[i];
                if( counters->data.i[j] > 1 )
                {
                    labels->data.i[i] = k;
                    counters->data.i[j]--;
                    counters->data.i[k]++;
                    break;
                }
            }

    __END__;

	cvConvert( centers, centers_arr );
    cvReleaseMat( &centers );
    cvReleaseMat( &old_centers );
    cvReleaseMat( &counters );
}


#if defined(WIN32)
void cvStartTimer(CvTimer * t){
	SYSTEMTIME st;
	FILETIME ft;
	GetSystemTime(&st);
	SystemTimeToFileTime(&st,&ft);
	t = ft.dwLowDateTime;
}
double cvTimerElapsed(CvTimer * t){
	SYSTEMTIME st;
	FILETIME ft;
	DWORD res;
	GetSystemTime(&st);
	SystemTimeToFileTime(&st,&ft);
	res = ft.dwLowDateTime - t;
	return ((double) res)/10000000.0;
}
#else
void cvStartTimer(CvTimer * t){
	gettimeofday(t, NULL);
}
double cvTimerElapsed(CvTimer * t){
	struct timeval now;
	double res; 
	gettimeofday(&now,NULL);
	timersub(&now, t, &now);
	res = now.tv_sec + now.tv_usec/1000000.0;
	return res;
}
#endif
