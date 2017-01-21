#ifndef __CV_EXT_H
#define __CV_EXT_H

#include <cxcore.h>
#include <cv.h>

//////  Additional Data Structures & functions on structures
#define CV_MAT_ROW_PTR( mat, row, type )   \
	((type *)((mat).data.ptr + (size_t)(mat).step*(row)))

typedef struct CvCircle32f {
	CvPoint2D32f c;
	float r;
} CvCircle32f;
CV_INLINE CvCircle32f cvCircle32f(float x, float y, float rad){
	CvCircle32f c;
	c.c = cvPoint2D32f(x,y);
	c.r = rad;
	return c;
}

typedef struct CvRect32f {
	float x;
	float y;
	float width;
	float height;
};

CV_INLINE CvRect32f cvRect32f(float x, float y, float width, float height){
	CvRect32f r;
	r.x = x;
	r.y = y;
	r.width = width;
	r.height = height;
	return r;
}

CV_INLINE CvFont cvFont( double scale CV_DEFAULT(1.0) ){
	CvFont font;
	cvInitFont(&font, CV_FONT_HERSHEY_PLAIN, scale, scale, 0, 1, CV_AA);
	return font;
}

/** Hash table implementation */
#include <glib.h>
typedef GHashTable CvHashTable;

/** Create hash table data structure.
 * @param key_type data type for hash table key 
 * @param value_type data type for hash table data
 */
CVAPI(CvHashTable *) cvCreateHashTable( CvTypeInfo * key_type, CvTypeInfo * value_type);

/** Returns number of elements in hash table */
CVAPI(int) cvHashTableSize(CvHashTable * h);

/** Insert an item into the hash table.  No copying is performed.
 * @param h     The hash table
 * @param key   Key used to find data.
 * @param value Data asociated with key
 */
CVAPI(void) cvHashTableInsert(CvHashTable * h, void * key, void * value);

/** Find data in hash table by its associated key */
CVAPI(void *) cvHashTableLookup(CvHashTable * h, const void * key);

CVAPI(void) cvReleaseHashTable( CvHashTable ** h );

/** Timer implementation */
#if defined(WIN32)
#include <olectl.h>
typedef DWORD CvTimer;
#else
#include <stdlib.h>
#include <sys/time.h>
typedef struct timeval CvTimer;
#endif 
void cvStartTimer(CvTimer * t);
double cvTimerElapsed(CvTimer * t);

/** 
 * Create image pyramid -- arguments are similar to cvCreateImage
 * @param size      Image width and height.
 * @param depth     Bit depth of image elements. 
 * @param nchannels Number of channels per element(pixel).
 * @param steps     Number of levels in the image pyramid
 * @ret An array of pointers to IplImages, such that the size(arr[i]) = size(arr[i-1])/2
 */
CVAPI(IplImage **) cvCreateImagePyr( CvSize size, int depth, int nchannels, int steps );

/**
 * Release image pyramid created with cvCreateImagePyr
 */
CVAPI(void) cvReleaseImagePyr( IplImage *** _pyr );

/**
 * calculate full image pyramid for a given image
 * @param src     Source image -- size mxn
 * @param dst     Destination pyramid -- must be base size m/2, n/2
 * @param filter  Type of the filter used for convolution;
 */
void cvPyrDownAll( CvArr * src, CvArr ** dst, int filter CV_DEFAULT(CV_GAUSSIAN_5x5) );
void cvPyrUpAll( CvArr * src, CvArr ** dst, int filter CV_DEFAULT(CV_GAUSSIAN_5x5) );


//////  Additional Functions ////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif

/** reallocate array(if neccessary) to match parameters
 * if parameters do not match, data is released and reallocated
 * if header is NULL, a new header is allocated and returned
 * @param header header to check against parameters 
 * @param rows   Number of rows in the matrix.
 * @param cols   Number of columns in the matrix.
 * @param type   Type of the matrix elements.
 * @see cvCreateMat
 */
CVAPI(CvMat *) cvReCreateMat( CvMat * header, int rows, int cols, int type );

// Return total number of elements ( rows x cols x nch ) for a given image or matrix
CVAPI(int) cvGetNumElements( CvArr * arr );

// Show scaled image -- computes minimum and maximum of image and scales between 0 and 255 before display
// If min and max are given, then these values are used for scaling
CVAPI(void) cvImageSC(const char * win, CvArr * arr, double min CV_DEFAULT(0), double max CV_DEFAULT(0));

// Show many images
// CVAPI(void) cvShowManyImages(char* title, int nArgs, ...);
CVAPI(void) cvShowManyImages(char* title, CvArr ** images, int nimages ); 

// Save array in Matlab 5 array format
CVAPI(int) cvSaveMat5(const char * fname, const CvArr * arr);


////////// File operations /////////////////////////////////

/** 
 * return only filename portion of path 
 * @param  name input path
 * @return      pointer into name indicating the filename portion of path
 */
CVAPI(char *) cvBaseName(const char * name);

/** 
 * extract filename without extension 
 * @param name  input path
 * @param dest  memory buffer for result, must be >= strlen(name)
 * @return      the filename with the .extension removed
 */
CVAPI(char *) cvRootName(const char * name, char * dest);

/** 
 * return only directory portion of path
 * @param name  input path
 * @param dest  memory buffer for result, must be >= strlen(name)
 * @return      the directory portion of the filename
 */
CVAPI(char *) cvDirName(const char * name, char * dest);



// Misc Image Operations /////////////////////////////////////////////////////////////////////////////

CVAPI(CvScalar) cvMedian( const CvArr * arr );
CVAPI(CvScalar) cvSelect( const CvArr * arr, int n );

/** reduces source array by performing 'func' on all elements of particular dimension,
 * thereby eliminating said dimension */
CV_EXTERN_C_FUNCPTR( CvScalar (*CvReduceFunc) (const CvArr *) );
CVAPI(void) cvReducePtr( const CvArr * arr, CvArr * dst, CvReduceFunc func, int dim CV_DEFAULT(0) );

/** 
 * Normalize pixels of image
 * The minimum and maximum of src are computed and src is scaled and shifted such that dest_min
 * and dest_max are the new bounds of the pixel values.  If mask!=NULL and mask(x,y)==0, then
 * (x,y) is not considered in the min/max computation and is not scaled in the destination image.
 * Inplace operation is also supported.
 * @param src       input array, any type, must be 1-channel
 * @param dst       output array, can be differet type from src unless mask!=NULL, 
 * @param dest_min  minimum of resulting array
 * @param dest_max  maximum of resulting array
 * @param mask      operation mask ( 8-bit 1 channel ).
 */
CV_INLINE void cvNormalizeImage(const CvArr * src, CvArr * dst, 
		double dest_min, double dest_max, const CvArr * mask CV_DEFAULT(NULL))
{
	cvNormalize( src, dst, dest_min, dest_max, CV_C, mask );
}

CVAPI(void)
cvKMeans3( const CvArr* samples_arr, int cluster_count,
           CvArr * centers_arr, CvArr* labels_arr, CvTermCriteria termcrit );
/** 
 * Compute the local minimum/maximum at a given point
 * @param arr    input array
 * @param p      point that defines local neighborhood
 * @param win    size of local neighborhood
 * @param min    returned minimum value in local neighborhood
 * @param max    returned maximum value in local neighborhood
 * @param minpos returned position of minimum value
 * @param maxpos returned position of maximum value
 */
void cvLocalMinMax(CvArr * arr, CvPoint p, CvSize win, double *min, double *max, CvPoint *minpos CV_DEFAULT(NULL), CvPoint *maxpos CV_DEFAULT(NULL));

/** Find the local min/max in arr, where (x,y) is a local minima if foreach x',y' in size x size
 * neighborhood around (x,y) arr(x,y) < arr(x', y') */
CVAPI(void) cvIsLocalMin( CvArr * arr, CvArr * dst, int size=3, CvArr * mask=NULL);
CVAPI(void) cvIsLocalMax( CvArr * arr, CvArr * dst, int size=3, CvArr * mask=NULL);

/** Find the points that are local min/maxima in arr, 
 *  (x,y) is a local minima if 
 *  mask(x,y) != 0 
 *  AND
 *  foreach x',y' in size x size  neighborhood around (x,y) 
 *  arr(x,y) < arr(x', y') AND arr(x,y) < arr_cmp(x', y') 
 *  @param arr          input array
 *  @param arr_cmp      array to compare against
 *  @param dst          CV_8UC1 output array
 *  @param test_center  whether or not to compare arr(x,y) to arr_cmp(x,y)
 *  @param size         width of 'local' neighborhood
 *  @param mask         operation mask
 **/
CVAPI(void) cvIsLocalMin2( CvArr * arr, CvArr * arr_cmp, CvArr * dst, int test_center=0, int size=3, CvArr * mask=NULL);
CVAPI(void) cvIsLocalMax2( CvArr * arr, CvArr * arr_cmp, CvArr * dst, int test_center=0, int size=3, CvArr * mask=NULL);

/* swap the elements of two arrays -- must be same size and nch */
#define SWAP_MAT_ELEMENTS(A, B, tmp)                         \
for(int _j=0; _j<A->rows; _j++){                             \
	for(int _i=0; _i<A->cols; _i++){                         \
		CV_SWAP(CV_MAT_ELEM( *A, typeof(tmp), _j, _i),       \
                CV_MAT_ELEM( *B, typeof(tmp), _j, _i), tmp); \
	}                                                        \
}

/** Determine if two arrays are identical 
 * @return 1 if a[i]=b[i] for all i 0 otherwise 
 */
CVAPI(int) cvIsEqual(CvArr * a_arr, CvArr * b_arr);

/** Return box that encapsulates the position, width, height, angle of image moments */
CvBox2D cvBox2DFromMoments( CvMoments * moments );

/** 
 * Computes the top N values in the array
 * @param arr input array
 * @param num N
 * @param min the minimum of the top N
 * @param max the maximum of the top N
 * @param minp location of the minimum value
 * @param maxp location of the maximum value
 * @param invert invert=1 indicates that calculate the minimum N values.
 */
CVAPI(void) cvMinMaxTopN( CvArr * arr, int num, float * min, float * max, 
		                  CvPoint * minp=NULL, CvPoint * maxp=NULL, int invert=0 );

/**
 * Rotate and zoom image around a point using cvGetQuadrangleSubPix
 * @param src    input image
 * @param dest   output image (can be different type)
 * @param center point to rotate/zoom about -- this becomes the center of the destination image
 * @param angle  angle to rotate the image by
 * @param scale  zoom factor
 */
CVAPI(void) cvRotateZoomImage(const CvArr * src, CvArr * dest, CvPoint2D32f center, float angle, 
		                      float scale);

/**
 * Rotate and zoom image around a point using cvWarpAffine
 * @param src     input image 
 * @param dest    output image (must be same type as input)
 * @param center  point to rotate/zoom about -- this becomes the center of the destination image
 * @param angle   angle to rotate the image by
 * @param scale   zoom factor
 * @param flags   see flags for cvWarpAffine
 * @param fillval value to use for points off the image
 */
CVAPI(void) cvRotateZoomImage2(const CvArr * src, CvArr * dest, CvPoint2D32f center, float angle, 
		                       float scale, int flags, CvScalar fillval CV_DEFAULT(cvScalarAll(0)));

/** compute gradient of image or array.
 * @param src
 * @param dx
 * @param dy
 * @param xwidth
 * @param ywidth
 */
CVAPI(void) cvGradient(const CvArr * src, CvArr * dst_x, CvArr * dst_y, int aperature_size CV_DEFAULT(3) );

/** compute difference of gaussian */
CV_INLINE void cvDoG(const CvArr * src, CvArr * dest, CvArr * temp, int param1, int param2){
	cvSmooth(src, temp, CV_GAUSSIAN, param1, param2, 0);
	cvSmooth(src, dest, CV_GAUSSIAN, param1+2, param2+2, 0);
	cvAbsDiff(dest,temp,dest);
}

// Misc Math Operations /////////////////////////////////////////////////////////////////////////////

/** Calculate angle between two points in space */
CV_INLINE float calcAngle(const CvPoint2D32f & p1, const CvPoint2D32f & p2){
	return atan2( p2.y-p1.y, p2.x-p1.x );
}

/** Calculate distance between two points in space */
CV_INLINE float L2(CvPoint2D32f p1, CvPoint2D32f p2){
	return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

// to fix problem with cvCeil on em64t and gcc 3.4
CV_INLINE  int  cvCeil1( double value )
{
	int temp = cvRound(value);
	Cv32suf diff;
	diff.f = (float)(temp - value);
	return temp + (diff.i < 0);
}

/**
 * Rearrange the quadrants of Fourier image so that the origin is at the image center 
 * @param src_arr 2-channel floating point array of fourier coefficients
 * @param dst_arr same size and shape as src_arr
 */
CVAPI(void) cvShiftDFT(CvArr * src_arr, CvArr * dst_arr );

/**	Returns the value ln[ (n-1)! ] for xx > 0. */
double cvLogGamma(double xx);

/** Returns the value ln( n! ) as a Floating-point number. */
double cvLogFactorial(int n);

/** 
 * Computes the binomial coefficient 
 * ( n )       n!
 * (   ) =  ---------
 * ( k )    (n-k)! k!
 */
long cvBinomial(int n, int k);

/**
 * Compute cross correlation coefficient between two arrays
 * Cross correlation of two arrays x,y is defined as 
 * SUM[ (x(i) - m_x) * (y(i)-m_y) ] / sqrt( SUM[ (x(i)-m_x)^2 ] * SUM[ (y(i)-m_y)^2 ] )
 * where
 * x(i) and y(i) are i'th elements of x and y arrays and
 * m_x and m_y are the means of x and y
 * @param src1 first input array
 * @param src2 second input array (same type as src1)
 * @param tmp1 temporary array (same type as src1)
 * @param tmp2 temporary array (same type as src2)
 * @return     cross correlation of src1 and src2
 */
CVAPI(double) cvCrossCorrelation(CvArr * src1, CvArr * src2, CvArr * tmp1, CvArr * tmp2);


// Linear Algebra functions /////////////////////////////////////////////////////////////////////

/** fill vec with given range of numbers */
CVAPI(CvArr *) cvRange(CvArr * vec, double start, double end);

/**
 * Compute the cholesky decomposition of a symmetric positive definite matrix A=U^T U
 * @param src input array
 * @param dst output array
 */
CVAPI(int) cvCholesky(CvArr * src, CvArr * dst);


/**
 * Project a vector into the subspace defined by basis
 * @param vec   input vector to project (Nx1)
 * @param basis MxN matrix containing the basis vectors as rows
 * @param proj  resulting projection
 */
CVAPI(void) cvProject( CvArr * vec, CvArr * basis, CvArr * proj );

/** 
 * Computes an orthonormal basis for the given sub-space
 * @param basis      MxN matrix of basis vectors as rows
 * @param orthoBasis orthonormal basis for given basis
 */
CVAPI(void) cvGramSchmidt( CvArr * basis, CvArr * orthoBasis );

#if 0
#define CV_PCA_DATA_AS_ROW 0 
#define CV_PCA_DATA_AS_COL 1
#define CV_PCA_USE_AVG 2
#define CV_PCA_RET_SDV 4
#define CV_PCA_RET_EVAL 0

/**
 * Perform Principal Components Analysis
 * @param data       CV_32F MxN input array with M input samples of N dimensions each
 * @param avg_mat    computed average of input rows -- must have N dimensions
 * @param eigenvals  computed eigenvalues associated with each eigenvector -- must have K dimensions
 *                   where K = MIN(M,N).  Can also be NULL
 * @param eigenvects computed eigenvectors -- each row is a vector -- must be KxN
 * @flags CV_PCA_DATA_AS_ROW indicates that each row in data is a sample
 *        CV_PCA_DATA_AS_COL indicates that each column in data is a sample
 *        CV_PCA_USE_AVG     does not compute the average of input data -- uses input avg_mat
 *        CV_PCA_RET_SDV     return standard deviations instead of eigen values
 *        CV_PCA_RET_EVAL    return eigen values of covariance matrix
 */
CVAPI(void) cvCalcPCA( CvArr * data, CvArr * avg_mat, CvArr * eigenvals, CvArr * eigenvects, 
		               int flags CV_DEFAULT( CV_PCA_DATA_AS_ROW | CV_PCA_RET_EVAL ) );

/** 
 * Project input data onto PCA subspace
 * @param data          CV_32F Nx1 input array
 * @param mean          CV_32F Nx1 mean computed in cvCalcPCA
 * @param eigenVectors  eigen vectors computed in cvCalcPCA
 * @param outResult     Kx1 projected data
 * @param dim           number eigen vectors to use in projection (-1 indicates ALL)
 */
CVAPI(void) cvProjectPCA( CvArr * data, CvArr * mean, CvArr * eigenVectors, CvArr * outResult, 
		                  int dim CV_DEFAULT(-1) );

/** 
 * 'Un'project projected data back to original space
 * @param eigenData  coefficients computed in cvCalcPCA
 * @param mean       mean compute in cvCalcPCA
 * @param eigenVectors eigenvectors computed in cvCalcPCA
 * @param outResult  unprojected data
 */
CVAPI(void) cvBackProjectPCA( CvArr * eigenData, CvArr * mean, CvArr * eigenVectors, CvArr * outResult );

#endif

// Miscellaneous OpenCV utility functions ///////////////////////////////////////////////////////
/** 
 * Error function for use with cvRedirectError to cause segfault instead of exit.  
 * Useful when debugging in gdb to get a backtrace rather than exiting.  
 * Call cvRedirectError( cvSegReport ) on the first line of your main() to enable this
 */
int cvSegReport( int status, const char * func_name,
		const char *err_msg, const char * file_name,
		int line, void *userdata );

#ifdef __cplusplus
} // extern "C"
#endif


#endif //__CV_EXT_H
