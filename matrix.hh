#ifndef __MATRIX_HH
#define __MATRIX_HH

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "cvmat.h"
#include "array_ops.hh"
//#include <gsl/gsl_linalg.h>
#include "vector3.h"

// from cvext.h
#ifdef __cplusplus
extern "C" {
#endif
	char * cvBaseName(const char * name);
#ifdef __cplusplus
}
#endif

// enable compilation of C++ Matrix functions
#ifndef USE_CPP_MATRIX
#define USE_CPP_MATRIX
#endif

// stupid template tricks to get correct type flag
template <typename T> struct MatrixDepth { static int depth() { assert(0); return 0; } };
template <> struct MatrixDepth<double> { static int depth() { return CV_64F; }};
template <> struct MatrixDepth<float>  { static int depth() { return CV_32F; }};
template <> struct MatrixDepth<int>    { static int depth() { return CV_32S; }};
template <> struct MatrixDepth<ushort> { static int depth() { return CV_16U; }};
template <> struct MatrixDepth<short>  { static int depth() { return CV_16S; }};
template <> struct MatrixDepth<uchar>  { static int depth() { return CV_8U; }};
template <> struct MatrixDepth<char>   { static int depth() { return CV_8S; }};
template <typename T, int nch> struct MatrixType { 
	static int type() { 
		return CV_MAKETYPE(MatrixDepth<T>::depth(), nch);
	}
};

#define STRINGIFY(a) #a
#define TOSTRING(a) STRINGIFY(a)
#define MATRIX_FROM_CVMAT( cvarr, matrix, func )\
do{\
	CvMat cvarr ## __cvmat;\
	cvGetMat(cvarr, & cvarr ## __cvmat);\
	switch( CV_MAT_TYPE(  cvarr ## __cvmat.type ) ){\
		case CV_8UC1:\
		{\
			Matrix<uchar,1> * matrix = (Matrix<uchar,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		case CV_8SC1:\
		{\
			Matrix<char,1> * matrix = (Matrix<char,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		case CV_16SC1:\
		{\
			Matrix<short,1> * matrix = (Matrix<short,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		case CV_32SC1:\
		{\
			Matrix<int,1> * matrix = (Matrix<int,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		case CV_32FC1:\
		{\
			Matrix<float,1> * matrix = (Matrix<float,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		case CV_64FC1:\
		{\
			Matrix<double,1> * matrix = (Matrix<double,1> *) & cvarr ## __cvmat;\
			func;\
		}\
		default:\
		{\
			char str[256];\
			sprintf(str, "CvMat type %d cannot be converted to Matrix", CV_MAT_TYPE(cvarr ## __cvmat.type));\
			cvError(CV_StsUnsupportedFormat, "MATRIX_FROM_CVMAT", str, __FILE__, __LINE__);\
		}\
	}\
}\
while( cvarr != cvarr )  // stupid compiler emits warning about while(0)

template <class T, int nch=1>
class Matrix : public CvMAT{
    
public:
/*    Matrix(int r, int c, bool alloc_data=true) : CvMAT(r,c,this->calcType,{
        cvInitMatHeader(this, r, c, this->calcType());
        if(alloc_data){
            cvCreateData(this);
        }
    }
*/
	Matrix(){
		type = 0;
	}
    Matrix(int r, int c) : CvMAT(r,c,calcType()) {
		cvZero(this);
	}
    Matrix(int r, int c, T * data) : CvMAT(r,c,this->calcType(), data) {}
    Matrix(const Matrix<T, nch> & m) : CvMAT((CvMAT) m){}
	Matrix(const CvMAT & m) : CvMAT(m) { /*assert(this->calcType())==*/}
	explicit Matrix( const _CvMAT_COPY_& mat_copy ) : CvMAT( mat_copy ) {}
    Matrix(const Vector3<T> & v) : CvMAT(3,1,this->calcType()){
		(*this) = v;
	}

	template <typename func_t>
	Matrix<T> operator=(Array_BIN_OP<Matrix<T>, func_t> op){
		op.func(*op.m_A, *op.m_B, (*this));
		return (*this);
	}  
	
	template <typename func_t>
	Matrix(Array_BIN_OP<Matrix<T>, func_t> op): CvMAT(op.m_A->width, op.m_A->height){
		op.func(*op.m_A, *op.m_B, (*this));
		//return (*this);
	}  
	// annoying redefinitions
    /*Matrix( const _CvMAT_T_& mat_t ) : CvMAT(mat_t) {}
    Matrix( const _CvMAT_INV_& inv_mat ) : CvMAT(inv_mat) {}
    Matrix( const _CvMAT_ADD_& mat_add ) : CvMAT(mat_add) {}
    Matrix( const _CvMAT_ADD_EX_& mat_add ) : CvMAT(mat_add) {}
    Matrix( const _CvMAT_SCALE_& scale_mat ) : CvMAT(scale_mat) {}
    Matrix( const _CvMAT_SCALE_SHIFT_& scale_shift_mat ) : CvMAT(scale_shift_mat) {}
    Matrix( const _CvMAT_MUL_& mmul ) : CvMAT(mmul) {}
    Matrix( const _CvMAT_MUL_ADD_& mmuladd ) : CvMAT(mmuladd) {}
    Matrix( const _CvMAT_LOGIC_& mat_logic ) : CvMAT(mat_logic) {}
    Matrix( const _CvMAT_UN_LOGIC_& mat_logic ) : CvMAT(mat_logic) {}
    Matrix( const _CvMAT_NOT_& not_mat ) : CvMAT(not_mat) {}
    Matrix( const _CvMAT_COPY_& mat_copy ) : CvMAT(mat_copy) {}
    Matrix( const _CvMAT_CVT_& mat_copy ) : CvMAT(mat_copy) {}
    Matrix( const _CvMAT_DOT_OP_& dot_mat ) : CvMAT(dot_mat) {}
    Matrix( const _CvMAT_SOLVE_& solve_mat ) : CvMAT(solve_mat) {}
    Matrix( const _CvMAT_CMP_& cmp_mat ) : CvMAT(cmp_mat) {}*/
	/*void init(int rows, int cols, T * data){
		cvInitMatHeader(this, rows, cols, calcType(), T);
	}*/
   	void operator= (const Matrix<T> & m){
		(*this) = (CvMat) m;
    }

	// conversion to Matrix class from IplImage
	static Matrix<T,nch> safe_cast( IplImage & im ){
		CvMat mat;
		cvGetMat(&im, &mat);
		return Matrix<T,nch>::safe_cast(mat); 
	}

	// checked casting to Matrix class from general CvMat
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
	// checked casting to Matrix class from general CvMat
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

	bool load(const char * fname, const char * name=NULL){
		return this->open(fname, name);
	}
	bool open(const char * fname, const char * name=NULL){
		const char * ext = fname+strlen(fname)-3;
		if(strcmp(ext, "xml")==0 || strcmp(ext, "yml")==0){
			return cv_open(fname, name);
		}
		if(strcmp(ext, "cvm")==0){
			return cvm_open(fname);
		}
		fprintf(stderr, "Matrix::open -- Extension '%s' not recognized\n", ext);
		return false;
	}

	bool cv_open(const char * fname, const char * name=NULL){
		CvFileStorage* fs = cvOpenFileStorage( fname, 0, CV_STORAGE_READ );
		CvMat * me = NULL; 
		CvFileNode * fn;
		if(name){
			fn = cvGetFileNodeByName(fs, NULL, name);
		}
		else{
			fn = cvGetFileNodeByName(fs,NULL, "CVMAT");
		}
		if(!fn){
			// find first object 
			fn = cvGetRootFileNode( fs );
			while(fn && CV_NODE_TYPE(fn->tag)!=CV_NODE_MAP){
				fn = (CvFileNode *) cvGetSeqElem(fn->data.seq, 0);
			}
		}
		if(fn){     
			CvAttrList attrlist;
			me = (Matrix<T,nch> *) cvRead( fs, fn, &attrlist);
			assert( CV_MAT_CN(me->type) == nch);
			assert( CV_MAT_DEPTH(me->type) == MatrixDepth<T>::depth() );
			if(me){
				//this->realloc(me->rows,me->cols);
				//cvCopy(me, this);
				(*this) = (*me);
				assert(this->rows == me->rows);
				assert(this->cols == me->cols);
				cvDecRefData(me);
				cvReleaseMat( &me );
				me = this;
			}
		}
		cvReleaseFileStorage( &fs );

		return me!=NULL;
	}

	/** load the image data from the given stream.
	 *  @param is an input stream that has already been opened and is ready to be read
	*/
	bool cvm_open(const char * fname){
        std::ifstream is(fname, std::ifstream::in | std::ifstream::binary);
        if(!is.is_open()) return false;
        bool ret = this->cvm_open(is);
        is.close();
        return ret;
	}
	bool cvm_open(std::istream & is){
		if(is.eof()) return false;
		
		CvMat header;
		unsigned long length;
		unsigned long offset;
		offset = is.tellg();
		//std::cout<<is.tellg()<<std::endl;
		is.seekg (0, std::ios::end);
		//std::cout<<is.tellg()<<std::endl;
		length = is.tellg();
		length-= offset;
		is.seekg (offset);
		//std::cout<<is.tellg()<<std::endl;	
		if(length<sizeof(CvMat)){
			std::cerr<<"Error reading CvMat Header -- expected size of"<<sizeof(CvMat)<<" got "<<length<<std::endl;
			return false;
		}
		
		// read header
		is.read((char *) &header, sizeof(CvMat));
		
		//std::cout<<is.tellg()<<std::endl;	

		if(length<sizeof(CvMat)+header.rows*header.step){
			std::cerr<<"Error reading array data -- not large enough -- expected "<<sizeof(CvMat)+header.rows*header.step<<" got "<<length<<std::endl;
			std::cerr<<"rows="<<header.rows<<" cols="<<header.cols<<" step="<<header.step<<std::endl;
			return false;
		}
			
		// blank out pointer fields
		header.data.ptr = 0;
		
		if(header.rows==0){
			std::cerr<<"Array has zero rows"<<std::endl;
			return false;
		}

		if(CV_MAT_TYPE(header.type) == this->calcType() ){
			int cstep;
			// allocate memory if neccessary
			this->realloc(header.rows, header.cols);

			// if the matrix is one row long, step is zero
			cstep = this->step;
			if(cstep==0){
				cstep = CV_ELEM_SIZE(this->type)*this->cols;
			}


			//printf("%d %d %d\n", this->rows, this->cols, this->step);
			
			//assert(this->step == header.step || fprintf(stderr, "this->step=%d header.step=%d dim=%d %d\n", this->step, header.step, header.rows, header.cols)<0);

			// read pixels
			is.read((char *)this->data.ptr, this->rows*cstep);
		}
		else{

			cvCreateData(&header);

			// if the matrix is one row long, step is zero
			//int cstep = this->step;
			int cstep = header.step;
			if(cstep==0){
				cstep = CV_ELEM_SIZE(this->type)*this->cols;
			}

			//printf("create(ipl_open): %p %d %d\n", header.imageData, header.width, header.height);
			is.read((char *)header.data.ptr, header.rows*header.step);
			
			// try to reshape to fit my number of channels
			if( (CV_MAT_CN(header.type)*header.rows*header.cols)%nch == 0 ){
				CvMat header2;
				if(header.rows % nch == 0){
					cvReshape(&header, &header2, nch, header.rows/nch);
				}
				else{
					cvReshape(&header, &header2, nch, header.rows);
				}
				this->realloc(header2.rows, header2.cols);
				cvScale(&header2, this);
			}
			else{
				this->realloc(header.rows, header.cols);
				// convert to my format
				cvScale(&header, this);
			}

			//printf("release(ipl_open): %p %d %d\n", this->imageData, this->width, this->height);
			cvReleaseData(&header);
		}
		//std::cout<<(*this)<<std::endl;
		return true;
	}
	bool cv_save(const char * fname) const{
		CvFileStorage* fs = cvOpenFileStorage( fname, 0, CV_STORAGE_WRITE );
		cvWrite( fs, "CVMAT", this, cvAttrList(0,0) );
		cvReleaseFileStorage( &fs );
		return true;
	}
	bool save(const char * fname) const{
        int len = strlen(fname);
        if(strcmp(fname+(len-3), "xml")==0 ||
           strcmp(fname+(len-3), "yml")==0) return cv_save(fname);
		if(strcmp(fname+(len-3), "cvm")==0) return cvm_save(fname);
        return false;
    }
    bool cvm_save(const char *fname) const{
        // the binary flag is neccessary on windows!
        std::ofstream os(fname, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
        if(!os.is_open()) return false;
        bool ret = this->cvm_save(os);
        os.close();
        return ret;
    }   
    bool cvm_save(std::ostream & os) const{
        //std::cout<<os.tellp()<<std::endl;
		int cstep = this->cols*CV_ELEM_SIZE(this->type);
		if(CV_IS_MAT_CONT( this->type )){
	        os.write((char *) this, sizeof(CvMat));
	        os.write((char *) this->data.ptr, this->rows*cstep);
		}
		// non-contiguous, must write in pieces and modify header
		else{
			CvMat header;
			cvInitMatHeader(&header, this->rows, this->cols, this->type);
			//printf("Matrix::cvm_save() -- header.step=%d this->step=%d\n", header.step, this->step);
			os.write((char *) &header, sizeof(CvMat));
			assert(header.step == cstep);
			for(int i=0; i<this->rows; i++){
				os.write((char *)(this->data.ptr+i*this->step), cstep);
			}
		}
        return true;
    }  
	/*Matrix<T> & operator = (const Vector3<T> & v){
	  assert(rows>=3);
		(*this)(0,0) = v[0];
		(*this)(1,0) = v[1];
		(*this)(2,0) = v[2];
		return (*this);
	}*/
    /* copying and filling with a constant */
    Matrix<T,nch> & operator = ( const CvMAT& mat ){
		CvMAT::operator=(mat);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const CvMat& mat ){
		CvMAT::operator=(mat);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const IplImage& img ){
		CvMAT::operator=(img);
		return (*this);
	}
    CvMAT & operator = ( double fillval ){
		CvMAT::operator=(fillval);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const CvScalar& fillval ){
		CvMAT::operator=(fillval);
		return (*this);
	}
	
    
    /* b = op(a1, a2,...) */
    Matrix<T,nch> & operator = ( const _CvMAT_T_& mat_t ){
		CvMAT::operator=(mat_t);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const _CvMAT_INV_& inv_mat ){
		CvMAT::operator=(inv_mat);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const _CvMAT_ADD_& mat_add ){
		CvMAT::operator=(mat_add);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const _CvMAT_ADD_EX_& mat_add ){
		CvMAT::operator=(mat_add);
		return (*this);
	}
    Matrix<T,nch> & operator = ( const _CvMAT_SCALE_& scale_mat ){
		CvMAT::operator=(scale_mat);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_SCALE_SHIFT_& scale_shift_mat ){
		CvMAT::operator=(scale_shift_mat);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_MUL_& mmul ){
		CvMAT::operator=(mmul);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_MUL_ADD_& mmuladd ){
		CvMAT::operator=(mmuladd);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_LOGIC_& mat_logic ){
		CvMAT::operator=(mat_logic);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_UN_LOGIC_& mat_logic ){
		CvMAT::operator=(mat_logic);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_NOT_& not_mat ){
		CvMAT::operator=(not_mat);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_DOT_OP_& dot_mul ){
		CvMAT::operator=(dot_mul);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_SOLVE_& solve_mat ){
		CvMAT::operator=(solve_mat);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_CMP_& cmp_mat ){
		CvMAT::operator=(cmp_mat);
		return (*this);
	}
    Matrix<T> & operator = ( const _CvMAT_CVT_& mat_cvt ){
		CvMAT::operator=(mat_cvt);
		return (*this);
	}
	Matrix<T> & operator = (const Vector3<T> & v){
		(*this)(0,0) = v[0];
		(*this)(1,0) = v[1];
		(*this)(2,0) = v[2];
		return (*this);
	}

    /* copy matrix data, not only matrix header */
    CvMAT& operator = ( const _CvMAT_COPY_& mat_copy ){
		CvMAT::operator=(mat_copy);
		return (*this);
	}

	T dot(const Matrix<T> & m){
		return cvDotProduct(this,&m);
	}
    void zero(){
        cvZero(this);
    }
    void identity(){
        cvSetIdentity(this, cvScalar(1));
    }
    int calcType(){
        return MatrixType<T,nch>::type();
    }
	
    /*CvMatrix3 getMatrix3(){
        CvMatrix3 m;
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                m.m[i][j] = (*this)(i,j);
            }
        }
        return m;
    }*/

	Vector3<T> & rowAsVector3(const int & row){
		return *(Vector3<T> *)(this->data.ptr+row*this->step);
	}
	Vector2<T> & rowAsVector2(const int & row){
		return *(Vector2<T> *)(this->data.ptr+row*this->step);
	}
    T & operator() (const int &row, const int &col){
        return ((T*)(this->data.ptr+row*this->step))[col*nch];
    }
    T & operator() (const int &row, const int &col, const int & ch){
        return ((T*)(this->data.ptr+row*this->step))[col*nch+ch];
    }
    const T & operator() (const int &row, const int &col, const int & ch) const{
        return ((T*)(this->data.ptr+row*this->step))[col*nch+ch];
    }
	const T & get(const int &row, const int &col) const{
        return ((T*)(this->data.ptr+row*this->step))[col*nch];
    }
    T * operator [] (const int &row){
		return ((T*)(this->data.ptr+row*this->step));
	}
    const T * operator [] (const int &row) const{
		return ((T*)(this->data.ptr+row*this->step));
	}
	Vector3<T> toVector3() const{
        assert(this->rows >=3);
        return Vector3<T> ((*this)(0,0), (*this)(1,0), (*this)(2,0));
    }
   	Matrix<T> operator() (const Matrix<int> & indices){
		Matrix<T> r(indices.rows, this->cols);
		for(int i=0; i<indices.rows; i++){
			T* srow = (*this)[ indices[i][0] ];
			T* drow = r[i];
			for(int j=0; j<this->cols; j++){
				drow[j] = srow[j];
			}
		}
		return r;
	}
    const T & operator() (const int &row, const int &col) const{
        #ifdef DEBUG
            assert(row<this->rows);
            assert(row>=0);
            assert(col<this->cols);
            assert(col>=0);
            assert(cvPtr2D(this,row,col)==this->data.ptr+row*this->step+col*sizeof(T));
        #endif
        return ((T*)(this->data.ptr+row*this->step))[col*nch];
    }
	uchar * getData() { return this->data.ptr; }
	const uchar * getData() const { return this->data.ptr; }
	
	/** resize will reallocate data if neccessary */
	void realloc(const int &r, const int &c){
		if(this->rows*this->cols == r*c){
			this->reshape(r,c);
		}
		else{
			cvDecRefData(this);
			cvInitMatHeader(this, r, c, this->calcType());
			cvCreateData(this);	
		}
	}
	
	/** reshape array to different # of rows/columns -- new rows*cols must == old rows*cols */
    void reshape(const int &r, const int &c){
		assert(r*c == this->rows*this->cols || fprintf(stderr, "Old=%dx%d New=%dx%d\n", this->rows, this->cols, r, c)<0);
		cvReshape(this, this, 0, r);
	}
	
	void reshape(const int &r) {
		cvReshape(this, this, 0, r);
	}
	Matrix<T,nch> asshape( int newr ) const {
		CvMAT m = (*this);
		m.reshape( nch, newr );
		return safe_cast( m );
	}
	template <int nch2>
	static Matrix<T,nch> asshape( const Matrix<T,nch2> &mat, int newr ){
		CvMAT m = mat;
		m.reshape( nch, newr );
		return safe_cast( m );
	}
	Matrix<T,nch> rowrange( int i, int j ) const{
		return safe_cast( CvMAT::rowrange(i,j) );
	}
	Matrix<T,nch> row( int i ) const{
		return safe_cast( CvMAT::row(i) );
	}
	Matrix<T,nch> colrange( int i, int j ) const{
		return safe_cast( CvMAT::colrange(i,j) );
	}
	Matrix<T,nch> col( int i ) const{
		return safe_cast( CvMAT::col(i) );
	}
	Matrix<T,nch> rect( const CvRect & r ) const{
		return safe_cast( CvMAT::rect(r) );
	}
	Matrix<T,nch> diag( int diag = 0 ) const{
		return safe_cast( CvMAT::diag( diag ) );
	}
	/*PointArray<T> operator* (PointArray<double> &a){
        PointArray<double> b(a.rows);
        cvGEMM(this, &a, 1.0, NULL, 0.0, &b, CV_GEMM_A_T + CV_GEMM_B_T);
        return b;
    }*/
    /*Vector3<T> operator* (const Vector3<double> &v) const{
        return Vector3<double> (
            v(0)*cvmGet(this, 0, 0) + v(1)*cvmGet(this, 1, 0) + v(2)*cvmGet(this, 2, 0),
            v(0)*cvmGet(this, 0, 1) + v(1)*cvmGet(this, 1, 1) + v(2)*cvmGet(this, 2, 1),
            v(0)*cvmGet(this, 0, 2) + v(1)*cvmGet(this, 1, 2) + v(2)*cvmGet(this, 2, 2));
    }*/
    /*static void repeat( const CvMat* srcarr, CvMat* dstarr )
    {

        int pix_size = 8;
        int x, y, k, l;

        for( y = 0, k = 0; y < dstarr->height; y++ )
        {
            for( x = 0; x < dstarr->width; x += srcarr->width )
            {
                l = srcarr->width;
                if( l > dstarr->width - x )
                        l = dstarr->width - x;
                memcpy( dstarr->data.ptr + y*dstarr->step + x*pix_size,
                        srcarr->data.ptr + k*srcarr->step, l*pix_size );
            }
            if( ++k == srcarr->height )
                    k = 0;
        }
    } */  
    static void setRotation(CvMat * m, int axis, float angle){
        double c = cos(angle);
        double s = sin(angle);
        cvSetIdentity(m);
        switch(axis){
        case 0: //X
            cvmSet(m,1,1,c);
            cvmSet(m,1,2,s);
            cvmSet(m,2,1,-s);
            cvmSet(m,2,2,c);
            break;
        case 1: //Y
            cvmSet(m,0,0,c);
            cvmSet(m,0,2,-s);
            cvmSet(m,2,0,s);
            cvmSet(m,2,2,c);
            break;
        case 2: //Z
            cvmSet(m,0,0,c);
            cvmSet(m,0,1,s);
            cvmSet(m,1,0,-s);
            cvmSet(m,1,1,c);
            break;
        }
    }

    void rotateObject(int axis, float angle){
        Matrix<T>::rotateObject(this, axis, angle);
    }
    void rotateAxes(int axis, float angle){
        Matrix<T>::rotateAxes(this, axis, angle);
    }
   	static void rotateObject(CvMat *m, int axis, float angle){
		CvMat *r = cvCreateMat(3,3, m->type);
		CvMat *clone = cvCloneMat(m);
		setRotation(r, axis, -angle);
		cvGEMM(r, clone, 1.0, 0, 0.0, m);
		cvReleaseMat(&clone);
		cvReleaseMat(&r);
	}
    static void rotateAxes(CvMat * m, int axis, float angle){
        CvMat *r = cvCreateMat(3,3, m->type);
        CvMat *clone = cvCloneMat(m);
        setRotation(r, axis, angle);
        cvGEMM(r, clone, 1.0, 0, 0.0, m);
        cvReleaseMat(&clone);
        cvReleaseMat(&r);
    }
    void set(T * arr){
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                cvmSet(this, i, j, arr[i*cols+j]);
            }
        }
    }
	T sum(){
		return (T) cvSum(this).val[0];
	}
    ~Matrix(){
		//~CvMAT();
        cvDecRefData( this );
    }
};

template <class T, int nch>
inline std::ostream & print_pixel(std::ostream & out, const Matrix<T,nch> &m, const int & i, const int & j){
	out<<"(";
	for(int k=0; k<nch-1; k++){
		out<<m(i,j,k)<<", ";
	}
	out<<m(i,j,nch-1)<<")";
	return out;
}
template <class T>
inline std::ostream & print_pixel(std::ostream & out, const Matrix<T,1> &m, const int & i, const int & j){
	out<<m(i,j);
	return out;
}

template <class T, int nch>
inline std::ostream & print_row(std::ostream & out, const Matrix<T,nch> &m, const int & i){
    out<<"[";
    for(int j=0; j<(m.cols-1); j++){
		print_pixel(out, m, i, j);
		out<<", ";
    }
	print_pixel(out, m, i, m.cols-1);
    out<<"]";
    return out;
}

template <class T, int nch>
inline std::ostream & operator <<(std::ostream & out, const Matrix<T,nch> &m) {
    out<<"[";
    for(int i=0; i<(m.rows-1); i++){
        print_row(out, m, i);
        out<<",\n";
    }
    print_row(out, m, m.rows-1);
    out<<"]";
    return out;
}

template <class T>
Vector3<T> operator* (const Matrix<T> & m, const Vector3<T> & v){
	Vector3<T> r = Vector3<T> (v(0)*m(0,0) + v(1)*m(0,1) + v(2)*m(0,2),
					   v(0)*m(1,0) + v(1)*m(1,1) + v(2)*m(1,2),
					   v(0)*m(2,0) + v(1)*m(2,1) + v(2)*m(2,2));
					   return r;
}

typedef Matrix<double> dMat;
typedef Matrix<float>  fMat;
typedef Matrix<int>    iMat;
typedef Matrix<char>   cMat;

#include <cv/cvext.h>

#endif //MATRIX_HH
