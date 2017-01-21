#ifndef __CV_MATRIX_HPP
#define __CV_MATRIX_HPP

/*M///////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004 Roman Stanchak
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without lmatitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
// permit persons to whom the Software is furnished to do so, subject to the following 
// conditions:
//
// The above copyright notice and this permission notice shall be included in all copies 
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
// IN THE SOFTWARE.
//
//M*/

#include <cxcore.h>
#include <stdio.h>

namespace cv {

// stupid template tricks to get correct depth flag
template <typename T> struct MatrixDepth { static int depth() { assert(0); return 0; }};
template <> struct MatrixDepth<double> {	static int depth() { return CV_64F; }};
template <> struct MatrixDepth<float> {	static int depth() { return CV_32F; }};
template <> struct MatrixDepth<int> {	static int depth() { return CV_32S; }};
template <> struct MatrixDepth<short> {	static int depth() { return CV_16S; }};
template <> struct MatrixDepth<unsigned short> {	static int depth() { return CV_16S; }};
template <> struct MatrixDepth<char> {	static int depth() { return CV_8S; }};
template <> struct MatrixDepth<unsigned char> {	static int depth() { return CV_8U; }};

// forward declaration for matrix operator overloading defined in <cv/array_ops.hpp>
template <typename pixel_t, int nch> struct Array_OP;

/// C++ wrapper for CvMat structure
template <typename T, int nch=1> class Matrix : public CvMat {
public:
	
	/** Default Constructor */	
    Matrix(){
		_setParameters();
    }
	
    /** Constructor.  
	 * @param w width of new matrix
	 * @param h height of new matrix
	 */
    Matrix(int m, int n){
		_setParameters();
		this->realloc(m, n);
    }

	/** Copy constructor.
	 *  does not copy underlying data, just increments refcount */
	Matrix(const Matrix<T,nch> & mat ){
		//assert(mat.type == CV_MAT_TYPE( CV_MAKETYPE( MatrixDepth<T>::depth(), nch )) );
		cvInitMatHeader( this, mat.rows, mat.cols, type, mat.data.ptr, mat.step );
		this->refcount = mat.refcount;
		cvIncRefData(this);
	}

	/** Reshaped matrix */
	explicit inline Matrix( const Matrix<T,nch>& mat, int nrows ){
		type = 0;
		cvReshape(&mat, this, 0, nrows);
		this->refcount = mat.refcount;
		cvIncRefData( this );
	}

	/** Extract sub rectangle */
    explicit inline Matrix( const Matrix<T,nch>& mat, CvRect rect ){ /* submatrix */
		type = 0;
		cvGetSubArr( &mat, this, rect );
		this->refcount = mat.refcount;
		cvIncRefData( this );
	}

	/* submatrix:
	 *      k == 0 - i-th row
	 *      k > 0 - i-th column
	 *      k < 0 - i-th diagonal */
	explicit inline Matrix( const Matrix<T,nch>& mat, int k, int i ){
		type = 0;
		if( k == 0 )
			cvGetRow( &mat, this, i );
		else if( k > 0 )
			cvGetCol( &mat, this, i );
		else
			cvGetDiag( &mat, this, i );
		this->refcount = mat.refcount;
		cvIncRefData( this );
	}

	/** = operator */
	Matrix<T, nch> & operator=( const Matrix<T, nch> & mat ){
		this->release();
		cvInitMatHeader( this, mat.rows, mat.cols, type, mat.data.ptr, mat.step );
		this->refcount = mat.refcount;
		cvIncRefData( this );
		return *this;
	}

	/** =operator for +,-,*,/ overloads */
	Matrix<T,nch> & operator = (const Array_OP<T,nch> & A){
		A.assign(this);
		return (*this);
	}


	/** checked casting to Matrix class from general CvMat */
	static Matrix<T,nch> & safe_cast( CvMat & m ){
		assert(CV_MAT_CN(m.type)==nch || fprintf(stderr, "mat.nch=%d this->nch=%d\n", CV_MAT_CN(m.type), nch)<0);
		assert(CV_MAT_DEPTH(m.type)==MatrixDepth<T>::depth() || fprintf(stderr, "mat.depth=%d this->depth=%d\n", CV_MAT_DEPTH(m.type), MatrixDepth<T>::depth())<0);
		return *(Matrix<T,nch> *) &m;
	}

	static const Matrix<T,nch> & safe_cast( const CvMat & m ) {
		assert(CV_MAT_CN(m.type)==nch || fprintf(stderr, "mat.nch=%d this->nch=%d\n", CV_MAT_CN(m.type), nch)<0);
		assert(CV_MAT_DEPTH(m.type)==MatrixDepth<T>::depth() || fprintf(stderr, "mat.depth=%d this->depth=%d\n", CV_MAT_DEPTH(m.type), MatrixDepth<T>::depth())<0);
		return *(Matrix<T,nch> *) &m;
	}

	/** checked casting to Matrix class from general CvMat */
	static Matrix<T,nch> * safe_cast( CvMat * m ){
		assert(CV_MAT_CN(m->type)==nch);
		assert(CV_MAT_DEPTH(m->type)==MatrixDepth<T>::depth());
		return (Matrix<T,nch> *) m;
	}

	static const Matrix<T,nch> * safe_cast( const CvMat * m ) {
		assert(CV_MAT_CN(m->type)==nch);
		assert(CV_MAT_DEPTH(m->type)==MatrixDepth<T>::depth());
		return (Matrix<T,nch> *) m;
	}

   	/** function to copy contents of destination */
	Matrix<T,nch> & clone( const Matrix<T,nch> & mat){
		this->realloc(mat.rows, mat.cols);
		cvCopy( &mat, this );
		return * this;
	}

	/** set internal data buffer to given pointer */
	/*void setData(T * data, int w=-1, int h=-1){
		if(w!=-1 && h!=-1){
			cvInitMatHeader(this, h, w, CV_MAKETYPE(MatrixDepth<T>::depth(), nch));
		}
		this->data.ptr = (uchar *) data;	

		// increment refcount
		cvIncRefData( this );
	}*/

	/** return a pointer to the data at row i */
	T* operator[](const int &i) {
		return (T*) (this->data.ptr+(i*this->step));
	}

	/** return a constant pointer to the data at row i */
	const T* operator[](const int &i) const {
		return (T*) (this->data.ptr+(i*this->step));
	}
	
	/** const individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    inline const T & operator () (const int & x, const int & y, const int &ch=0) const {
		return  ((T*) (this->data.ptr+(y*this->step)))[x*nch+ch];
	}

	/** individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    inline T & operator () (const int & x, const int & y, const int &ch=0) {
		return  ((T*) (this->data.ptr+(y*this->step)))[x*nch+ch];
	}

	/** allocate or reallocate (if neccessary) the matrix buffer.  
	    No memory allocation will take place if w==matrix.width and h==matrix.height
		@param w width of the reallocd matrix
		@param h height of the reallocd matrix
	*/
    bool realloc(const int & m, const int & n){
		if(m<=0 || n<=0) return false;
        if(this->data.ptr){
            if(this->rows==m && 
               this->cols==n ){ 
                return false;
            }
			//printf("release(realloc): %p %d %d\n", this->data.ptr, this->width, this->height);
            //cvReleaseData(this);
			this->release();
        }
		_setParameters();
        cvInitMatHeader(this, m, n, CV_MAKETYPE(MatrixDepth<T>::depth(), nch));
        cvCreateData(this);
		//printf("create(realloc): %p %d %d\n", this->data.ptr, this->width, this->height);
		return true;
    }

	/** return the width of the matrix */
	const int & getWidth() const { return this->cols; }
	
	/** return the height of the matrix */
	const int & getHeight() const { return this->rows; }
	
	/** return the depth of the matrix */
	int getDepth() const { return CV_MAT_DEPTH(this->type); }
	
	/** return the number of channels in the matrix */
	int getNumChannels() const { return CV_MAT_CN(this->type); }

	/** row, column, sub-rectangle access */
	inline Matrix<T,nch> row( int r ) const
	{   return Matrix<T,nch>( *this, 0, r );  }


	inline Matrix<T,nch> col( int c ) const
	{   return Matrix<T,nch>( *this, 1, c );  }


	inline Matrix<T,nch> diag( int d ) const
	{   return Matrix<T,nch>( *this, -1, d );  }


	inline Matrix<T,nch> rect( CvRect rect ) const
	{   return Matrix<T,nch>( *this, rect );    }

	inline Matrix<T,nch> rowrange( int row1, int row2 ) const
	{   
		assert( 0 <= row1 && row1 < row2 && row2 <= height );
		return Matrix<T,nch>( *this, cvRect( 0, row1, width, row2 - row1 ));
	}

	inline Matrix<T,nch> colrange( int col1, int col2 ) const
	{   
		assert( 0 <= col1 && col1 < col2 && col2 <= width );
		return Matrix<T,nch>( *this, cvRect( col1, 0, col2 - col1, height ));
	}
	
	/** return header around same data with different dimensions */
	inline Matrix<T,nch> reshape( int rows ) const
	{
		return Matrix<T,nch>( *this, rows );
	}

    /* decrease matrix data reference counter and clear data pointer */
	inline void release()
	{   
		cvDecRefData( this );
	}

	bool load(const char * fname, const char * name=NULL){
		CvFileStorage* fs = cvOpenFileStorage( fname, 0, CV_STORAGE_READ );
		CvMat * me = NULL; 
		CvFileNode * fn;
		if(name){
			fn = cvGetFileNodeByName(fs, NULL, name);
		}
		else{
			//char str[16];
			//sscanf(cvBaseName(fname), "%[^.].%*s", str);
			//sprintf(str, "CVMAT");
			fn = cvGetFileNodeByName(fs,NULL, "CVMAT");
		}
		if(fn){     
			CvAttrList attrlist;
			me = (CvMat *) cvRead( fs, fn, &attrlist);
			if(me){
				assert( CV_MAT_CN(me->type) == nch);
				assert( CV_MAT_DEPTH(me->type) == MatrixDepth<T>::depth() );
				(*this) = (Matrix<T,nch>::safe_cast(*me));
				cvReleaseMat(&me);
				me = this;
			}
		}
		cvReleaseFileStorage( &fs );

		return me!=NULL;
	}
	bool save(const char * fname){
		//char str[16];
		//char * basename = cvBaseName( fname );
		//assert( sscanf(basename, "%[^.].%*s", str)!=0 );
		CvFileStorage* fs = cvOpenFileStorage( fname, 0, CV_STORAGE_WRITE );
		cvWrite( fs, "CVMAT", this, cvAttrList(0,0) );
		cvReleaseFileStorage( &fs );
		return true;
	}

	~Matrix(){
		this->release();
	}

private:
	void _setParameters(){
		// T mytypevar;
		// cvInitMatHeader(this, 0, 0, CV_MAKETYPE(MatrixDepth<T>::depth(), nch) );
		this->data.ptr=NULL;
		this->refcount=NULL;
	}
};

} // namespace cv

#include <cv/array_ops.tcc>
ARRAY_OPS_IMPL(cv::Matrix)

#endif //__CV_MATRIX_HPP
