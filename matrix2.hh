/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#ifndef _CVMAT_HPP_
#define _CVMAT_HPP_


#include <iostream>
#include <string.h>
#include <stdio.h>
#include <cxcore.h>

#undef min
#undef max

/****************************************************************************************\
*                                   C++ Matrix Class                                     *
\****************************************************************************************/
class Matrix : public CvMat
{
protected:

public:
    /* helper methods for retrieving/setting matrix elements */
    static double get( const uchar* ptr, int type, int coi = 0 );
    static void set( uchar* ptr, int type, int coi, double d );
    static void set( uchar* ptr, int type, int coi, int i );
    static void set( uchar* ptr, int type, double d );
    static void set( uchar* ptr, int type, int i );

    /******************* constructors ********************/
    /* empty */
    explicit Matrix();

    /* creation */
    explicit Matrix( int rows, int cols, int type, void* data, int step = CV_AUTOSTEP );
    explicit Matrix( int rows, int type, void* data, int step = CV_AUTOSTEP );
    explicit Matrix( int rows, int cols, int type );
    explicit Matrix( int rows, int type );
    
    /* extracting part of an existing matrix */
    explicit Matrix( const CvMat& mat, CvRect rect ); /* submatrix */
    explicit Matrix( const CvMat& mat, int k, int i ); /* submatrix:
                                                    k == 0 - i-th row
                                                    k > 0 - i-th column
                                                    k < 0 - i-th diagonal */
    /* copying */
    Matrix( const CvMat& mat );
    Matrix( const Matrix& mat );
    Matrix( const IplImage& img );

    /* Matrix b = op(a1,a2,...) */
	explicit Matrix( const _Matrix_T_& mat_t );
    explicit Matrix( const _Matrix_INV_& inv_mat );
    explicit Matrix( const _Matrix_ADD_& mat_add );
    explicit Matrix( const _Matrix_ADD_EX_& mat_add );
    explicit Matrix( const _Matrix_SCALE_& scale_mat );
    explicit Matrix( const _Matrix_SCALE_SHIFT_& scale_shift_mat );
    explicit Matrix( const _Matrix_MUL_& mmul );
    explicit Matrix( const _Matrix_MUL_ADD_& mmuladd );
    explicit Matrix( const _Matrix_LOGIC_& mat_logic );
    explicit Matrix( const _Matrix_UN_LOGIC_& mat_logic );
    explicit Matrix( const _Matrix_NOT_& not_mat );
    explicit Matrix( const _Matrix_COPY_& mat_copy );
    explicit Matrix( const _Matrix_CVT_& mat_copy );
    explicit Matrix( const _Matrix_DOT_OP_& dot_mul );
    explicit Matrix( const _Matrix_SOLVE_& solve_mat );
    explicit Matrix( const _Matrix_CMP_& cmp_mat );

    /* desctructor */
    ~Matrix();

    /* copying and filling with a constant */
    Matrix& operator = ( const Matrix& mat );
    Matrix& operator = ( const CvMat& mat );
    Matrix& operator = ( const IplImage& img );
    Matrix& operator = ( double fillval );
    Matrix& operator = ( const CvScalar& fillval );
    
    /* b = op(a1, a2,...) */
    Matrix& operator = ( const _Matrix_T_& mat_t );
    Matrix& operator = ( const _Matrix_INV_& inv_mat );
    Matrix& operator = ( const _Matrix_ADD_& mat_add );
    Matrix& operator = ( const _Matrix_ADD_EX_& mat_add );
    Matrix& operator = ( const _Matrix_SCALE_& scale_mat );
    Matrix& operator = ( const _Matrix_SCALE_SHIFT_& scale_shift_mat );
    Matrix& operator = ( const _Matrix_MUL_& mmul );
    Matrix& operator = ( const _Matrix_MUL_ADD_& mmuladd );
    Matrix& operator = ( const _Matrix_LOGIC_& mat_logic );
    Matrix& operator = ( const _Matrix_UN_LOGIC_& mat_logic );
    Matrix& operator = ( const _Matrix_NOT_& not_mat );
    Matrix& operator = ( const _Matrix_DOT_OP_& dot_mul );
    Matrix& operator = ( const _Matrix_SOLVE_& solve_mat );
    Matrix& operator = ( const _Matrix_CMP_& cmp_mat );
    Matrix& operator = ( const _Matrix_CVT_& mat_cvt );

    /* copy matrix data, not only matrix header */
    Matrix& operator = ( const _Matrix_COPY_& mat_copy );

    /* augmented assignments */
    Matrix& operator += ( const CvMat& mat );
    Matrix& operator += ( double val );
    Matrix& operator += ( const CvScalar& val );
    Matrix& operator += ( const _Matrix_SCALE_& scale_mat );
    Matrix& operator += ( const _Matrix_SCALE_SHIFT_& scale_mat );
    Matrix& operator += ( const _Matrix_MUL_& mmul );

    Matrix& operator -= ( const CvMat& mat );
    Matrix& operator -= ( double val );
    Matrix& operator -= ( const CvScalar& val );
    Matrix& operator -= ( const _Matrix_SCALE_& scale_mat );
    Matrix& operator -= ( const _Matrix_SCALE_SHIFT_& scale_mat );
    Matrix& operator -= ( const _Matrix_MUL_& mmul );

    Matrix& operator *= ( const CvMat& mat );
    Matrix& operator *= ( double val );
    Matrix& operator *= ( const CvScalar& val );
    Matrix& operator *= ( const _Matrix_SCALE_& scale_mat );
    Matrix& operator *= ( const _Matrix_SCALE_SHIFT_& scale_mat );

    Matrix& operator &= ( const CvMat& mat );
    Matrix& operator &= ( double val );
    Matrix& operator &= ( const CvScalar& val );

    Matrix& operator |= ( const CvMat& mat );
    Matrix& operator |= ( double val );
    Matrix& operator |= ( const CvScalar& val );

    Matrix& operator ^= ( const CvMat& mat );
    Matrix& operator ^= ( double val );
    Matrix& operator ^= ( const CvScalar& val );

    /* various scalar charactertics */
    double norm( int norm_type = CV_L2 ) const;
    double norm( CvMat& mat, int norm_type = CV_L2 ) const;
    CvScalar sum() const;

    double det() const;
//    double trace() const;

    _Matrix_T_  t() const; /* transposition */
    _Matrix_INV_ inv(int method = 0) const;
    /* inversion using one of the following methods:
          method = 0 - Gaussian elimination,
          method = 1 - SVD */

    _Matrix_DOT_OP_  mul( const Matrix& mat ) const;
    _Matrix_DOT_OP_  mul( const _Matrix_SCALE_& mat ) const;

    _Matrix_DOT_OP_  div( const Matrix& mat ) const;
    _Matrix_DOT_OP_  div( const _Matrix_SCALE_& mat ) const;

    _Matrix_DOT_OP_  min( const Matrix& mat ) const;
    _Matrix_DOT_OP_  max( const Matrix& mat ) const;
    _Matrix_DOT_OP_  min( double value ) const;
    _Matrix_DOT_OP_  max( double value ) const;
    double          min( CvPoint* minloc = 0 ) const;
    double          max( CvPoint* maxloc = 0 ) const;

    _Matrix_DOT_OP_  abs() const;
    
    /* accessing matrix elements */
    _MatrixElem_ operator ()( int row );
    _MatrixConstElem_ operator ()( int row ) const;

    _MatrixElem_ operator ()( int row, int col );
    _MatrixConstElem_ operator ()( int row, int col ) const;

    _MatrixElem_ operator ()( CvPoint loc );
    _MatrixConstElem_ operator ()( CvPoint loc ) const;

    _MatrixElemCn_ operator()( int row, int col, int coi );
    double operator()( int row, int col, int coi ) const;

    _MatrixElemCn_ operator()( CvPoint pt, int coi );
    double operator()( CvPoint pt, int coi ) const;

    void* ptr( int row );
    const void* ptr( int row ) const;

    void* ptr( int row, int col );
    const void* ptr( int row, int col ) const;

    void* ptr( CvPoint pt );
    const void* ptr( CvPoint pt ) const;

    /* accessing matrix parts */
    Matrix row( int row ) const;
    Matrix rowrange( int row1, int row2 ) const;
    Matrix col( int col ) const;
    Matrix colrange( int col1, int col2 ) const;
    Matrix rect( CvRect rect ) const;
    Matrix diag( int diag = 0 ) const;

    _Matrix_COPY_ clone() const;

    /* convert matrix */
    _Matrix_CVT_ cvt( int newdepth = -1, double scale = 1,
                     double shift = 0 ) const;

    /* matrix transformation */
    void  reshape( int newcn, int newrows = 0 );
    void  flipX();
    void  flipY();
    void  flipXY();

    /* matrix I/O: use dynamically linked runtime libraries */
    void  write( const char* name = 0, FILE* f = 0, const char* fmt = 0 );
    void  read( char** name = 0, FILE* f = 0 );

    /* decrease matrix data reference counter and clear data pointer */
    void  release();

    /* internal use only! */
    bool  istemp() const;
	friend std::ostream & operator<<(std::ostream & out, const Matrix & m);
	friend std::ostream & print_row(std::ostream & out, const Matrix & m, const int & i);

protected:
    
    void  create( int rows, int cols, int type );
};


/* !!! Internal Use Only !!! */
/* proxies for matrix elements */

/* const_A(i,j) */
struct   _MatrixConstElem_
{
    explicit _MatrixConstElem_( const uchar* ptr, int type );
    operator CvScalar () const;
    double operator ()( int coi = 0 ) const;

    uchar* ptr;
    int type;
};


/* A(i,j,cn) or A(i,j)(cn) */
struct   _MatrixElemCn_
{
    explicit _MatrixElemCn_( uchar* ptr, int type, int coi );
    operator double() const;
    
    _MatrixElemCn_& operator = ( const _MatrixConstElem_& elem );
    _MatrixElemCn_& operator = ( const _MatrixElemCn_& elem );
    _MatrixElemCn_& operator = ( const CvScalar& scalar );
    _MatrixElemCn_& operator = ( double d );
    _MatrixElemCn_& operator = ( float f );
    _MatrixElemCn_& operator = ( int i );

    uchar* ptr;
    int type;
};


/* A(i,j) */
struct   _MatrixElem_ : public _MatrixConstElem_
{
    explicit _MatrixElem_( uchar* ptr, int type );
    _MatrixElemCn_ operator ()( int coi = 0 );

    _MatrixElem_& operator = ( const _MatrixConstElem_& elem );
    _MatrixElem_& operator = ( const _MatrixElem_& elem );
    _MatrixElem_& operator = ( const _MatrixElemCn_& elem );
    _MatrixElem_& operator = ( const CvScalar& val );
    _MatrixElem_& operator = ( double d );
    _MatrixElem_& operator = ( float f );
    _MatrixElem_& operator = ( int i );
};


/* !!! Internal Use Only !!! */
/* temporary classes for matrix operations */
struct   _Matrix_TMP_ : public Matrix
{
    _Matrix_TMP_();
    _Matrix_TMP_( const _Matrix_TMP_& mat );
    _Matrix_TMP_& operator = ( const _Matrix_TMP_& mat );

    explicit _Matrix_TMP_( const _Matrix_T_& mat_t );
    explicit _Matrix_TMP_( const _Matrix_INV_& inv_mat );
    explicit _Matrix_TMP_( const _Matrix_ADD_& mat_add );
    explicit _Matrix_TMP_( const _Matrix_ADD_EX_& mat_add );
    explicit _Matrix_TMP_( const _Matrix_SCALE_& scale_mat );
    explicit _Matrix_TMP_( const _Matrix_SCALE_SHIFT_& scale_shift_mat );
    explicit _Matrix_TMP_( const _Matrix_MUL_& mmul );
    explicit _Matrix_TMP_( const _Matrix_MUL_ADD_& mmuladd );
    explicit _Matrix_TMP_( const _Matrix_LOGIC_& mmuladd );
    explicit _Matrix_TMP_( const _Matrix_UN_LOGIC_& mmuladd );
    explicit _Matrix_TMP_( const _Matrix_NOT_& not_mat );
    explicit _Matrix_TMP_( const _Matrix_DOT_OP_& dot_mul );
    explicit _Matrix_TMP_( const _Matrix_SOLVE_& solve_mat );
    explicit _Matrix_TMP_( const _Matrix_CMP_& cmp_mat );
    explicit _Matrix_TMP_( const _Matrix_CVT_& cmp_mat );
    
    /* extracting parts from Matrix */
    explicit _Matrix_TMP_( const _Matrix_TMP_& mat, CvRect rect );
    explicit _Matrix_TMP_( const _Matrix_TMP_& mat, int k, int i );
};


struct   _Matrix_BASE_OP_
{
    _Matrix_BASE_OP_() {};
    virtual operator _Matrix_TMP_() const = 0;

    _Matrix_DOT_OP_  mul( const Matrix& a ) const;
    _Matrix_DOT_OP_  mul( const _Matrix_SCALE_& a ) const;

    _Matrix_DOT_OP_  div( const Matrix& a ) const;
    _Matrix_DOT_OP_  div( const _Matrix_SCALE_& a ) const;

    _Matrix_DOT_OP_  max( const Matrix& a ) const;
    _Matrix_DOT_OP_  min( const Matrix& a ) const;

    _Matrix_DOT_OP_  max( double value ) const;
    _Matrix_DOT_OP_  min( double value ) const;

    double          max( CvPoint* maxloc = 0 ) const;
    double          min( CvPoint* minloc = 0 ) const;

    _Matrix_DOT_OP_  abs() const;

    _Matrix_INV_     inv( int method = 0 ) const;
    _Matrix_T_       t() const;

    _Matrix_TMP_     row( int row ) const;
    _Matrix_TMP_     rowrange( int row1, int row2 ) const;
    _Matrix_TMP_     col( int col ) const;
    _Matrix_TMP_     colrange( int col1, int col2 ) const;
    _Matrix_TMP_     rect( CvRect rect ) const;
    _Matrix_TMP_     diag( int diag = 0 ) const;
    _Matrix_CVT_     cvt( int newdepth = -1, double scale = 1, double shift = 0 ) const;
    
    double          norm( int norm_type = CV_L2 ) const;
    double          det() const;
 //   double          trace() const;
    CvScalar        sum() const;
	virtual ~_Matrix_BASE_OP_() {}
};


/* (A^t)*alpha */
struct   _Matrix_T_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_T_( const Matrix* a );
    explicit _Matrix_T_( const Matrix* a, double alpha );
    
    double det() const;
    double norm( int normType = CV_L2 ) const;
    operator _Matrix_TMP_() const;

    Matrix  a;
    double alpha;
};


/* inv(A) */
struct   _Matrix_INV_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_INV_( const Matrix* mat, int method );
    operator _Matrix_TMP_() const;

    Matrix a;
    int method;
};


/* (A^ta)*(B^tb)*alpha */
struct   _Matrix_MUL_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_MUL_( const Matrix* a, const Matrix* b, int t_ab );
    explicit _Matrix_MUL_( const Matrix* a, const Matrix* b,
                          double alpha, int t_abc );
    operator _Matrix_TMP_() const;

    Matrix* a;
    Matrix* b;
    double alpha;
    int t_ab; /* (t_ab & 1) = ta, (t_ab & 2) = tb */
};


/* (A^ta)*(B^tb)*alpha + (C^tc)*beta */
struct   _Matrix_MUL_ADD_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_MUL_ADD_( const Matrix* a, const Matrix* b,
                              const Matrix* c, int t_abc );
    explicit _Matrix_MUL_ADD_( const Matrix* a, const Matrix* b, double alpha,
                              const Matrix* c, double beta, int t_abc );
    operator _Matrix_TMP_() const;

    Matrix* a;
    Matrix* b;
    double alpha;
    Matrix* c;
	double beta;
    int t_abc; /* (t_abc & 1) = ta, (t_abc & 2) = tb, (t_abc & 4) = tc */
};


/* A + B*beta */
struct   _Matrix_ADD_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_ADD_( const Matrix* a, const Matrix* b, double beta = 1 );
    operator _Matrix_TMP_() const;

    double   norm( int norm_type = CV_L2 ) const;
    _Matrix_DOT_OP_ abs() const;

    Matrix* a;
    Matrix* b;
    double beta;
};


/* A*alpha + B*beta + gamma */
struct   _Matrix_ADD_EX_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_ADD_EX_( const Matrix* a, double alpha,
                             const Matrix* b, double beta, double gamma = 0 );
    operator _Matrix_TMP_() const;

    Matrix* a;
	double alpha;
    Matrix* b;
    double beta, gamma;
};


/* A*alpha */
struct   _Matrix_SCALE_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_SCALE_( const Matrix* a, double alpha );
    operator _Matrix_TMP_() const;

    _Matrix_DOT_OP_  mul( const Matrix& a ) const;
    _Matrix_DOT_OP_  mul( const _Matrix_SCALE_& a ) const;

    _Matrix_DOT_OP_  div( const Matrix& a ) const;
    _Matrix_DOT_OP_  div( const _Matrix_SCALE_& a ) const;

    Matrix* a;
    double alpha;
};


/* A*alpha + beta */
struct   _Matrix_SCALE_SHIFT_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_SCALE_SHIFT_( const Matrix* a, double alpha, double beta );
    operator _Matrix_TMP_() const;

    _Matrix_DOT_OP_  abs() const;

    Matrix* a;
    double alpha, beta;
};


/* (A & B), (A | B) or (A ^ B) */
struct   _Matrix_LOGIC_ : public _Matrix_BASE_OP_
{
    enum Op { AND = 0, OR = 1, XOR = 2 };
    explicit _Matrix_LOGIC_( const Matrix* a, const Matrix* b, Op op, int flags = 0 );
    operator _Matrix_TMP_() const;

    Matrix* a;
    Matrix* b;
    Op op;
    int flags;
};


/* (A & scalar), (A | scalar) or (A ^ scalar) */
struct   _Matrix_UN_LOGIC_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_UN_LOGIC_( const Matrix* a, double alpha,
                               _Matrix_LOGIC_::Op op, int flags = 0 );
    operator _Matrix_TMP_() const;

    Matrix* a;
    double alpha;
    _Matrix_LOGIC_::Op op;
    int flags;
};


/* ~A */
struct   _Matrix_NOT_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_NOT_( const Matrix* a );
    operator _Matrix_TMP_() const;

    Matrix* a;
};


/* conversion of data type */
struct   _Matrix_CVT_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_CVT_( const Matrix* a, int newdepth = -1,
                          double scale = 1, double shift = 0 );
    operator _Matrix_TMP_() const;

    Matrix a;
    int newdepth;
    double scale, shift;
};


/* conversion of data type */
struct   _Matrix_COPY_
{
    explicit _Matrix_COPY_( const Matrix* a );
    operator Matrix() const;
    Matrix* a;
};


/* a.op(b), where op = mul, div, min, max ... */
struct   _Matrix_DOT_OP_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_DOT_OP_( const Matrix* a, const Matrix* b,
                             int op, double alpha = 1 );
    operator _Matrix_TMP_() const;

    Matrix a; /* keep the left operand copy */
    Matrix* b;
    int op;
    double alpha;
};


/* A.inv()*B or A.pinv()*B */
struct   _Matrix_SOLVE_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_SOLVE_( const Matrix* a, const Matrix* b, int method );
    operator _Matrix_TMP_() const;

    Matrix* a;
    Matrix* b;
    int method;
};


/* A <=> B */
struct   _Matrix_CMP_ : public _Matrix_BASE_OP_
{
    explicit _Matrix_CMP_( const Matrix* a, const Matrix* b, int cmp_op );
    explicit _Matrix_CMP_( const Matrix* a, double alpha, int cmp_op );
    operator _Matrix_TMP_() const;

    Matrix* a;
    Matrix* b;
    double alpha;
    int cmp_op;
};


/************************* _MatrixConstElem_ inline methods ******************************/

inline _MatrixConstElem_::_MatrixConstElem_(const uchar* p, int t) : ptr((uchar*)p), type(t)
{}


inline _MatrixConstElem_::operator CvScalar() const
{
    CvScalar scalar;
    cvRawDataToScalar( ptr, type, &scalar );

    return scalar;
}


inline double _MatrixConstElem_::operator ()( int coi ) const
{   return Matrix::get( ptr, type, coi );    }


inline _MatrixElemCn_::_MatrixElemCn_( uchar* p, int t, int coi ) :
    ptr(p), type(CV_MAT_DEPTH(t))
{
    if( coi )
    {
        assert( (unsigned)coi < (unsigned)CV_MAT_CN(t) );
        ptr += coi * CV_ELEM_SIZE(type);
    }
}


inline _MatrixElemCn_::operator double() const
{   return Matrix::get( ptr, type ); }


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( const _MatrixConstElem_& elem )
{
    if( type == elem.type )
        memcpy( ptr, elem.ptr, CV_ELEM_SIZE(type) );
    else
    {
        assert( CV_MAT_CN(elem.type) == 1 );
        Matrix::set( ptr, type, 0, elem(0));
    }

    return *this;
}


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( const _MatrixElemCn_& elem )
{
    if( type == elem.type )
        memcpy( ptr, elem.ptr, CV_ELEM_SIZE(type) );
    else
        Matrix::set( ptr, type, 0, (double)elem );
    return *this;
}


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( const CvScalar& scalar )
{   
    Matrix::set( ptr, type, 0, scalar.val[0] );
    return *this;
}


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( double d )
{   
    Matrix::set( ptr, type, 0, d );
    return *this;
}


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( float f )
{   
    Matrix::set( ptr, type, 0, (double)f );
    return *this;
}


inline _MatrixElemCn_& _MatrixElemCn_::operator = ( int i )
{   
    Matrix::set( ptr, type, 0, i );
    return *this;
}


inline _MatrixElem_::_MatrixElem_( uchar* p, int t ) : _MatrixConstElem_( (const uchar*)p, t )
{}


inline _MatrixElemCn_ _MatrixElem_::operator ()( int coi )
{   return _MatrixElemCn_( ptr, type, coi ); }


inline _MatrixElem_& _MatrixElem_::operator = ( const _MatrixConstElem_& elem )
{
    if( type == elem.type )
        memcpy( ptr, elem.ptr, CV_ELEM_SIZE(type) );
    else
    {
        assert( CV_MAT_CN( type ^ elem.type ) == 0 );
        CvScalar sc = (CvScalar)elem;
        cvScalarToRawData( &sc, ptr, type, 0 );
    }

    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( const _MatrixElem_& elem )
{
    *this = (const _MatrixConstElem_&)elem;
    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( const _MatrixElemCn_& elem )
{
    if( type == elem.type )
        memcpy( ptr, elem.ptr, CV_ELEM_SIZE(type) );
    else
        Matrix::set( ptr, type, (double)elem );

    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( const CvScalar& scalar )
{
    cvScalarToRawData( &scalar, ptr, type, 0 );
    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( double d )
{
    Matrix::set( ptr, type, d );
    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( float f )
{
    Matrix::set( ptr, type, (double)f );
    return *this;
}


inline _MatrixElem_& _MatrixElem_::operator = ( int i )
{
    Matrix::set( ptr, type, i );
    return *this;
}


/********************************** Matrix inline methods ********************************/

inline Matrix::Matrix()
{
    memset( this, 0, sizeof(*this));
}


inline Matrix::Matrix( int rows, int cols, int type, void* data, int step )
{
    cvInitMatHeader( this, rows, cols, type, data, step );
}


inline Matrix::Matrix( int rows, int type, void* data, int step )
{
    cvInitMatHeader( this, rows, 1, type, data, step );
}


inline void Matrix::create( int rows, int cols, int type )
{
    cvInitMatHeader( this, rows, cols, type );
    cvCreateData( this );
	cvZero(this);
}


inline Matrix::Matrix( int rows, int cols, int type )
{
    create( rows, cols, type );
}


inline Matrix::Matrix( int rows, int type )
{
    create( rows, 1, type );
}


inline Matrix::Matrix( const CvMat& mat )
{
    memcpy( this, &mat, sizeof(mat));
    if( refcount )
        (*refcount)++;
}


inline Matrix::Matrix( const Matrix& mat )
{
    memcpy( this, &mat, sizeof(mat));
    if( refcount )
        (*refcount)++;
}


inline Matrix::Matrix( const IplImage& img )
{
    cvGetMat( &img, this );
}


inline void Matrix::release()
{
    if( type != 0 && refcount && --refcount[0] == 0 )
    {
        uchar* ptr = (uchar*)refcount + sizeof(refcount)*2;
        cvFree( (void**)&ptr );
    }
    data.ptr = 0;
    refcount = 0;
}


inline Matrix::~Matrix()
{
    release();
}


inline bool Matrix::istemp() const
{
    return (type & CV_MAT_TEMP_FLAG) != 0;
}


inline Matrix& Matrix::operator = ( const Matrix& mat )
{
    if( this != &mat )
    {
        release();
        memcpy( this, &mat, sizeof(mat));
        if( refcount )
            (*refcount)++;
    }
    return *this;
}


inline Matrix& Matrix::operator = ( const CvMat& mat )
{
    *this = (const Matrix&)mat;
    return *this;
}


inline Matrix& Matrix::operator = ( const IplImage& img )
{
    release();
    cvGetMat( &img, this );

    return *this;
}


inline Matrix& Matrix::operator = ( double fillval )
{
    cvSet( this, cvScalarAll(fillval) );
    return *this;
}


inline Matrix& Matrix::operator = ( const CvScalar& fillval )
{
    cvSet( this, fillval );
    return *this;
}


inline Matrix& Matrix::operator += ( const CvMat& mat )
{
    cvAdd( this, &mat, this );
    return *this;
}


inline Matrix& Matrix::operator += ( double val )
{
    cvAddS( this, cvScalar(val), this );
    return *this;
}


inline Matrix& Matrix::operator += ( const CvScalar& val )
{
    cvAddS( this, val, this );
    return *this;
}


inline Matrix& Matrix::operator -= ( const CvMat& mat )
{
    cvSub( this, &mat, this );
    return *this;
}


inline Matrix& Matrix::operator -= ( double val )
{
    cvSubS( this, cvScalar(val), this );
    return *this;
}


inline Matrix& Matrix::operator -= ( const CvScalar& val )
{
    cvSubS( this, val, this );
    return *this;
}


inline Matrix& Matrix::operator *= ( const CvMat& mat )
{
    cvMul( this, &mat, this );
    return *this;    
}


inline Matrix& Matrix::operator *= ( double val )
{
    cvScale( this, this, val, 0 );
    return *this;
}


inline Matrix& Matrix::operator *= ( const CvScalar& val )
{
    cvScaleAdd( this, val, 0, this );
    return *this;
}


inline Matrix& Matrix::operator &= ( const CvMat& mat )
{
    cvAnd( this, &mat, this );
    return *this;
}


inline Matrix& Matrix::operator &= ( double val )
{
    cvAndS( this, cvScalarAll(val), this );
    return *this;
}


inline Matrix& Matrix::operator &= ( const CvScalar& val )
{
    cvAndS( this, val, this );
    return *this;
}


inline Matrix& Matrix::operator |= ( const CvMat& mat )
{
    cvOr( this, &mat, this );
    return *this;
}


inline Matrix& Matrix::operator |= ( double val )
{
    cvOrS( this, cvScalarAll(val), this );
    return *this;
}


inline Matrix& Matrix::operator |= ( const CvScalar& val )
{
    cvOrS( this, val, this );
    return *this;
}


inline Matrix& Matrix::operator ^= ( const CvMat& mat )
{
    cvXor( this, &mat, this );
    return *this;
}


inline Matrix& Matrix::operator ^= ( double val )
{
    cvXorS( this, cvScalarAll(val), this );
    return *this;
}


inline Matrix& Matrix::operator ^= ( const CvScalar& val )
{
    cvXorS( this, val, this );
    return *this;
}


inline double Matrix::norm( int normType ) const
{   return cvNorm( this, 0, normType ); }


inline double Matrix::min( CvPoint* minloc ) const
{   
    double t;
    cvMinMaxLoc( this, &t, 0, minloc, 0, 0 );
    return t;
}

inline double Matrix::max( CvPoint* maxloc ) const
{   
    double t;
    cvMinMaxLoc( this, 0, &t, 0, maxloc, 0 );
    return t;
}


inline double Matrix::norm( CvMat& mat, int normType ) const
{   return cvNorm( this, &mat, normType );  }


inline CvScalar Matrix::sum() const
{   return cvSum( this );   }


inline double Matrix::det() const
{   return cvDet( this );   }


inline void Matrix::reshape( int newcn, int newrows )
{   cvReshape( this, this, newcn, newrows );    }


inline void Matrix::flipX()
{   cvFlip( this, this, 1 );    }


inline void Matrix::flipY()
{   cvFlip( this, this, 0 );    }


inline void Matrix::flipXY()
{   cvFlip( this, this, -1 );   }


inline _MatrixElem_ Matrix::operator ()( int row )
{   return _MatrixElem_( CV_MAT_ELEM_PTR( *this, row, 0 ), type );   }


inline _MatrixConstElem_ Matrix::operator ()( int row ) const
{   return _MatrixConstElem_( CV_MAT_ELEM_PTR( *this, row, 0 ), type ); }


inline _MatrixElem_ Matrix::operator ()( int row, int col )
{   return _MatrixElem_( CV_MAT_ELEM_PTR( *this, row, col ), type ); }


inline _MatrixConstElem_ Matrix::operator ()( int row, int col ) const
{   return _MatrixConstElem_( CV_MAT_ELEM_PTR( *this, row, col ), type ); }


inline _MatrixElemCn_ Matrix::operator()( int row, int col, int coi )
{   return _MatrixElemCn_( CV_MAT_ELEM_PTR( *this, row, col ), type, coi );  }


inline _MatrixElemCn_ Matrix::operator()( CvPoint pt, int coi )
{   return _MatrixElemCn_( CV_MAT_ELEM_PTR( *this, pt.y, pt.x ), type, coi );  }


inline double Matrix::operator()( int row, int col, int coi ) const
{   return get( CV_MAT_ELEM_PTR( *this, row, col ), type, coi );    }


inline _MatrixElem_ Matrix::operator ()( CvPoint pt )
{   return _MatrixElem_( CV_MAT_ELEM_PTR( *this, pt.y, pt.x ), type ); }


inline _MatrixConstElem_ Matrix::operator ()( CvPoint pt ) const
{   return _MatrixConstElem_( CV_MAT_ELEM_PTR( *this, pt.y, pt.x ), type ); }


inline double Matrix::operator()( CvPoint pt, int coi ) const
{   return get( CV_MAT_ELEM_PTR( *this, pt.y, pt.x ), type, coi );    }


inline void* Matrix::ptr( int row )
{   return CV_MAT_ELEM_PTR( *this, row, 0 );    }


inline const void* Matrix::ptr( int row ) const
{   return (const void*)CV_MAT_ELEM_PTR( *this, row, 0 );   }


inline void* Matrix::ptr( int row, int col )
{   return CV_MAT_ELEM_PTR( *this, row, col );  }


inline const void* Matrix::ptr( int row, int col ) const
{   return (const void*)CV_MAT_ELEM_PTR( *this, row, col ); }


inline void* Matrix::ptr( CvPoint pt )
{   return CV_MAT_ELEM_PTR( *this, pt.y, pt.x );    }


inline const void* Matrix::ptr( CvPoint pt ) const
{   return (const void*)CV_MAT_ELEM_PTR( *this, pt.y, pt.x ); }


inline _Matrix_INV_ Matrix::inv( int method ) const
{   return _Matrix_INV_( this, method ); }


inline _Matrix_T_ Matrix::t() const
{   return _Matrix_T_( this );   }


inline _Matrix_COPY_ Matrix::clone() const
{   return _Matrix_COPY_( this ); }

inline _Matrix_CVT_ Matrix::cvt( int newdepth, double scale, double shift ) const
{   return _Matrix_CVT_( this, newdepth, scale, shift ); }

inline Matrix::Matrix( const CvMat& mat, CvRect rect )
{   
	type = 0;
    cvGetSubArr( &mat, this, rect );
    cvIncRefData( this );
}


/* submatrix:
     k == 0 - i-th row
     k > 0 - i-th column
     k < 0 - i-th diagonal */
inline Matrix::Matrix( const CvMat& mat, int k, int i )
{
    type = 0;
    if( k == 0 )
        cvGetRow( &mat, this, i );
    else if( k > 0 )
        cvGetCol( &mat, this, i );
    else
        cvGetDiag( &mat, this, i );
    cvIncRefData( this );
}


inline Matrix Matrix::row( int r ) const
{   return Matrix( *this, 0, r );  }


inline Matrix Matrix::col( int c ) const
{   return Matrix( *this, 1, c );  }


inline Matrix Matrix::diag( int d ) const
{   return Matrix( *this, -1, d );  }


inline Matrix Matrix::rect( CvRect rect ) const
{   return Matrix( *this, rect );    }

inline Matrix Matrix::rowrange( int row1, int row2 ) const
{   
    assert( 0 <= row1 && row1 < row2 && row2 <= height );
    return Matrix( *this, cvRect( 0, row1, width, row2 - row1 ));
}

inline Matrix Matrix::colrange( int col1, int col2 ) const
{   
    assert( 0 <= col1 && col1 < col2 && col2 <= width );
    return Matrix( *this, cvRect( col1, 0, col2 - col1, height ));
}

inline _Matrix_DOT_OP_ Matrix::mul( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( this, &mat, '*' );   }

inline _Matrix_DOT_OP_ Matrix::mul( const _Matrix_SCALE_& mat ) const
{   return _Matrix_DOT_OP_( this, mat.a, '*', mat.alpha );   }

inline _Matrix_DOT_OP_ Matrix::div( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( this, &mat, '/' );  }

inline _Matrix_DOT_OP_ Matrix::div( const _Matrix_SCALE_& mat ) const
{   return _Matrix_DOT_OP_( this, mat.a, '/', 1./mat.alpha );    }

inline _Matrix_DOT_OP_ Matrix::min( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( this, &mat, 'm' );   }

inline _Matrix_DOT_OP_ Matrix::max( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( this, &mat, 'M' );   }

inline _Matrix_DOT_OP_ Matrix::min( double value ) const
{   return _Matrix_DOT_OP_( this, 0, 'm', value );   }

inline _Matrix_DOT_OP_ Matrix::max( double value ) const
{   return _Matrix_DOT_OP_( this, 0, 'M', value );   }

inline _Matrix_DOT_OP_ Matrix::abs() const
{   return _Matrix_DOT_OP_( this, 0, 'a', 0 );   }

/****************************************************************************************\
*                               binary operations (+,-,*)                                *
\****************************************************************************************/

/*
* PART I. Scaling, shifting, addition/subtraction operations
*/

/* (mat2^t) = (mat1^t) * scalar */
inline _Matrix_T_ operator * ( const _Matrix_T_& a, double alpha )
{   return _Matrix_T_( &a.a, a.alpha*alpha );  }

/* (mat2^t) = scalar * (mat1^t) */
inline _Matrix_T_ operator * ( double alpha, const _Matrix_T_& a )
{   return _Matrix_T_( &a.a, a.alpha*alpha );  }

/* -(mat^t) */
inline _Matrix_T_ operator - ( const _Matrix_T_& a )
{   return _Matrix_T_( &a.a, -a.alpha ); }

/* mat_scaled = mat * scalar */
inline _Matrix_SCALE_ operator * ( const Matrix& a, double alpha )
{   return _Matrix_SCALE_( &a, alpha );  }

/* mat_scaled = scalar * mat */
inline _Matrix_SCALE_ operator * ( double alpha, const Matrix& a )
{   return _Matrix_SCALE_( &a, alpha );  }

/* mat_scaled2 = mat_scaled1 * scalar */
inline _Matrix_SCALE_ operator * ( const _Matrix_SCALE_& a, double alpha )
{   return _Matrix_SCALE_( a.a, a.alpha*alpha ); }

/* mat_scaled2 = scalar * mat_scaled1 */
inline _Matrix_SCALE_ operator * ( double alpha, const _Matrix_SCALE_& a )
{   return _Matrix_SCALE_( a.a, a.alpha*alpha ); }

/* -mat_scaled */
inline _Matrix_SCALE_ operator - ( const _Matrix_SCALE_& a )
{   return _Matrix_SCALE_( a.a, -a.alpha ); }


/* mat_scaled_shifted = mat + scalar */
inline _Matrix_SCALE_SHIFT_ operator + ( const Matrix& a, double beta )
{   return _Matrix_SCALE_SHIFT_( &a, 1, beta );  }

/* mat_scaled_shifted = scalar + mat */
inline _Matrix_SCALE_SHIFT_ operator + ( double beta, const Matrix& a )
{   return _Matrix_SCALE_SHIFT_( &a, 1, beta );  }

/* mat_scaled_shifted = mat - scalar */
inline _Matrix_SCALE_SHIFT_ operator - ( const Matrix& a, double beta )
{   return _Matrix_SCALE_SHIFT_( &a, 1, -beta ); }

/* mat_scaled_shifted = scalar - mat */
inline _Matrix_SCALE_SHIFT_ operator - ( double beta, const Matrix& a )
{   return _Matrix_SCALE_SHIFT_( &a, -1, beta ); }

/* mat_scaled_shifted = mat_scaled + scalar */
inline _Matrix_SCALE_SHIFT_ operator + ( const _Matrix_SCALE_& a, double beta )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, beta );   }

/* mat_scaled_shifted = scalar + mat_scaled */
inline _Matrix_SCALE_SHIFT_ operator + ( double beta, const _Matrix_SCALE_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, beta );   }

/* mat_scaled_shifted = mat_scaled - scalar */
inline _Matrix_SCALE_SHIFT_ operator - ( const _Matrix_SCALE_& a, double beta )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, -beta );  }

/* mat_scaled_shifted = scalar - mat_scaled */
inline _Matrix_SCALE_SHIFT_ operator - ( double beta, const _Matrix_SCALE_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, -a.alpha, beta );  }

/* mat_scaled_shifted2 = mat_scaled_shifted1 + scalar */
inline _Matrix_SCALE_SHIFT_ operator + ( const _Matrix_SCALE_SHIFT_& a, double beta )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, a.beta + beta );  }

/* mat_scaled_shifted2 = scalar + mat_scaled_shifted1 */
inline _Matrix_SCALE_SHIFT_ operator + ( double beta, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, a.beta + beta );  }

/* mat_scaled_shifted2 = mat_scaled_shifted1 - scalar */
inline _Matrix_SCALE_SHIFT_ operator - ( const _Matrix_SCALE_SHIFT_& a, double beta )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha, a.beta - beta );  }

/* mat_scaled_shifted2 = scalar - mat_scaled_shifted1 */
inline _Matrix_SCALE_SHIFT_ operator - ( double beta, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, -a.alpha, beta - a.beta ); }

/* mat_scaled_shifted2 = mat_scaled_shifted1 * scalar */
inline _Matrix_SCALE_SHIFT_ operator * ( const _Matrix_SCALE_SHIFT_& a, double alpha )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha*alpha, a.beta*alpha ); }

/* mat_scaled_shifted2 = scalar * mat_scaled_shifted1 */
inline _Matrix_SCALE_SHIFT_ operator * ( double alpha, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, a.alpha*alpha, a.beta*alpha ); }

/* -mat_scaled_shifted */
inline _Matrix_SCALE_SHIFT_ operator - ( const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_SCALE_SHIFT_( a.a, -a.alpha, -a.beta ); }


/* -mat1 */
inline _Matrix_SCALE_ operator - ( const Matrix& a )
{   return _Matrix_SCALE_( &a, -1 );   }

/* mat_add = mat1 + mat2 */
inline _Matrix_ADD_ operator + ( const Matrix& a, const Matrix& b )
{   return _Matrix_ADD_( &a, &b );   }

/* mat_add = mat1 - mat2 */
inline _Matrix_ADD_ operator - ( const Matrix& a, const Matrix& b )
{   return _Matrix_ADD_( &a, &b, -1 );    }

/* mat_add = mat_scaled1 + mat2 */
inline _Matrix_ADD_ operator + ( const _Matrix_SCALE_& a, const Matrix& b )
{   return _Matrix_ADD_( &b, a.a, a.alpha );  }

/* mat_add = mat1 + mat_scaled2 */
inline _Matrix_ADD_ operator + ( const Matrix& b, const _Matrix_SCALE_& a )
{   return _Matrix_ADD_( &b, a.a, a.alpha );  }

/* -mat_add */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_ADD_& a )
{   return _Matrix_ADD_EX_( a.a, -1, a.b, -a.beta ); }

/* mat_add = mat_scaled1 - mat2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_& a, const Matrix& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, &b, -1 ); }

/* mat_add = mat1 - mat_scaled2 */
inline _Matrix_ADD_ operator - ( const Matrix& b, const _Matrix_SCALE_& a )
{   return _Matrix_ADD_( &b, a.a, -a.alpha ); }

/* mat_add = mat_scaled_shifted1 + mat2 */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_SCALE_SHIFT_& a, const Matrix& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, &b, 1, a.beta ); }

/* mat_add = mat1 + mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator + ( const Matrix& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_ADD_EX_( a.a, a.alpha, &b, 1, a.beta ); }

/* mat_add = mat_scaled_shifted1 - mat2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_SHIFT_& a, const Matrix& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, &b, -1, a.beta ); }

/* mat_add = mat1 - mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator - ( const Matrix& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_ADD_EX_( a.a, -a.alpha, &b, 1, -a.beta ); }

/* mat_add = mat_scaled_shifted1 + mat_scaled2 */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_SCALE_SHIFT_& a, const _Matrix_SCALE_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, b.alpha, a.beta ); }

/* mat_add = mat_scaled1 + mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_SCALE_& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, b.alpha, a.beta ); }

/* mat_add = mat_scaled_shifted1 - mat_scaled2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_SHIFT_& a, const _Matrix_SCALE_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, -b.alpha, a.beta ); }

/* mat_add = mat_scaled1 - mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_ADD_EX_( a.a, -a.alpha, b.a, b.alpha, -a.beta ); }

/* mat_add = mat_scaled1 + mat_scaled2 */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_SCALE_& a, const _Matrix_SCALE_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, b.alpha ); }

/* mat_add = mat_scaled1 - mat_scaled2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_& a, const _Matrix_SCALE_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, -b.alpha ); }

/* mat_add = mat_scaled_shifted1 + mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_SCALE_SHIFT_& a,
                                const _Matrix_SCALE_SHIFT_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, b.alpha, a.beta + b.beta ); }

/* mat_add = mat_scaled_shifted1 - mat_scaled_shifted2 */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_SCALE_SHIFT_& a,
                                const _Matrix_SCALE_SHIFT_& b )
{   return _Matrix_ADD_EX_( a.a, a.alpha, b.a, -b.alpha, a.beta - b.beta ); }

/* mat_add2 = mat_add1 + scalar */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_ADD_EX_& a, double gamma )
{   return _Matrix_ADD_EX_( a.a, a.alpha, a.b, a.beta, a.gamma + gamma ); }

/* mat_add2 = scalar + mat_add1 */
inline _Matrix_ADD_EX_ operator + ( double gamma, const _Matrix_ADD_EX_& a )
{   return _Matrix_ADD_EX_( a.a, a.alpha, a.b, a.beta, a.gamma + gamma ); }

/* mat_add2 = mat_add1 - scalar */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_ADD_EX_& a, double gamma )
{   return _Matrix_ADD_EX_( a.a, a.alpha, a.b, a.beta, a.gamma - gamma ); }

/* mat_add2 = scalar - mat_add1 */
inline _Matrix_ADD_EX_ operator - ( double gamma, const _Matrix_ADD_EX_& a )
{   return _Matrix_ADD_EX_( a.a, -a.alpha, a.b, -a.beta, gamma - a.gamma ); }

/* mat_add2 = mat_add1 * scalar */
inline _Matrix_ADD_EX_ operator * ( const _Matrix_ADD_EX_& a, double alpha )
{   return _Matrix_ADD_EX_( a.a, a.alpha*alpha, a.b, a.beta*alpha, a.gamma*alpha ); }

/* mat_add2 = scalar * mat_add1 */
inline _Matrix_ADD_EX_ operator * ( double alpha, const _Matrix_ADD_EX_& a )
{   return _Matrix_ADD_EX_( a.a, a.alpha*alpha, a.b, a.beta*alpha, a.gamma*alpha ); }

/* mat_add2 = mat_add1 + scalar */
inline _Matrix_ADD_EX_ operator + ( const _Matrix_ADD_& a, double gamma )
{   return _Matrix_ADD_EX_( a.a, 1, a.b, a.beta, gamma ); }

/* mat_add2 = scalar + mat_add1 */
inline _Matrix_ADD_EX_ operator + ( double gamma, const _Matrix_ADD_& a )
{   return _Matrix_ADD_EX_( a.a, 1, a.b, a.beta, gamma ); }

/* mat_add2 = mat_add1 - scalar */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_ADD_& a, double gamma )
{   return _Matrix_ADD_EX_( a.a, 1, a.b, a.beta, -gamma ); }

/* mat_add2 = scalar - mat_add1 */
inline _Matrix_ADD_EX_ operator - ( double gamma, const _Matrix_ADD_& a )
{   return _Matrix_ADD_EX_( a.a, -1, a.b, -a.beta, gamma ); }

/* mat_add2 = mat_add1 * scalar */
inline _Matrix_ADD_EX_ operator * ( const _Matrix_ADD_& a, double alpha )
{   return _Matrix_ADD_EX_( a.a, alpha, a.b, a.beta*alpha, 0 ); }

/* mat_add2 = scalar * mat_add1 */
inline _Matrix_ADD_EX_ operator * ( double alpha, const _Matrix_ADD_& a )
{   return _Matrix_ADD_EX_( a.a, alpha, a.b, a.beta*alpha, 0 ); }

/* -mat_add_ex */
inline _Matrix_ADD_EX_ operator - ( const _Matrix_ADD_EX_& a )
{   return _Matrix_ADD_EX_( a.a, -a.alpha, a.b, -a.beta, -a.gamma ); }


/*
* PART II. Matrix multiplication.
*/

/* mmul = mat1 * mat2 */
inline _Matrix_MUL_ operator * ( const Matrix& a, const Matrix& b )
{   return _Matrix_MUL_( &a, &b, 0 );    }

/* mmul = (mat1^t) * mat2 */
inline _Matrix_MUL_ operator * ( const _Matrix_T_& a, const Matrix& b )
{   return _Matrix_MUL_( &a.a, &b, a.alpha, 1 );   }

/* mmul = mat1 * (mat2^t) */
inline _Matrix_MUL_ operator * ( const Matrix& b, const _Matrix_T_& a )
{   return _Matrix_MUL_( &b, &a.a, a.alpha, 2 );  }

/* mmul = (mat1^t) * (mat2^t) */
inline _Matrix_MUL_ operator * ( const _Matrix_T_& a, const _Matrix_T_& b )
{   return _Matrix_MUL_( &a.a, &b.a, a.alpha*b.alpha, 3 );  }

/* mmul = mat_scaled1 * mat2 */
inline _Matrix_MUL_ operator * ( const _Matrix_SCALE_& a, const Matrix& b )
{   return _Matrix_MUL_( a.a, &b, a.alpha, 0 ); }

/* mmul = mat1 * mat_scaled2 */
inline _Matrix_MUL_ operator * ( const Matrix& b, const _Matrix_SCALE_& a )
{   return _Matrix_MUL_( &b, a.a, a.alpha, 0 ); }

/* mmul = (mat1^t) * mat_scaled1 */
inline _Matrix_MUL_ operator * ( const _Matrix_T_& a, const _Matrix_SCALE_& b )
{   return _Matrix_MUL_( &a.a, b.a, a.alpha*b.alpha, 1 ); }

/* mmul = mat_scaled1 * (mat2^t) */
inline _Matrix_MUL_ operator * ( const _Matrix_SCALE_& b, const _Matrix_T_& a )
{   return _Matrix_MUL_( b.a, &a.a, a.alpha*b.alpha, 2 ); }

/* mmul = mat_scaled1 * mat_scaled2 */
inline _Matrix_MUL_ operator * ( const _Matrix_SCALE_& a, const _Matrix_SCALE_& b )
{   return _Matrix_MUL_( a.a, b.a, a.alpha*b.alpha, 0 ); }

/* mmul2 = mmul1 * scalar */
inline _Matrix_MUL_ operator * ( const _Matrix_MUL_& a, double alpha )
{   return _Matrix_MUL_( a.a, a.b, a.alpha*alpha, a.t_ab ); }

/* mmul2 = scalar * mmul1 */
inline _Matrix_MUL_ operator * ( double alpha, const _Matrix_MUL_& a )
{   return _Matrix_MUL_( a.a, a.b, a.alpha*alpha, a.t_ab ); }

/* -mmul */
inline _Matrix_MUL_ operator - ( const _Matrix_MUL_& a )
{   return _Matrix_MUL_( a.a, a.b, -a.alpha, a.t_ab ); }

/* mmuladd = mmul + mat */
inline _Matrix_MUL_ADD_ operator + ( const _Matrix_MUL_& a, const Matrix& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b, 1, a.t_ab ); }

/* !!! Comment this off because of ambigous conversion error !!!
   mmuladd = mat + mmul */
/* inline _Matrix_MUL_ADD_ operator + ( const Matrix& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b, 1, a.t_ab ); }*/

/* mmuladd = mmul - mat */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_MUL_& a, const Matrix& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b, -1, a.t_ab ); }

/* !!! Comment this off because of ambigous conversion error !!!
   mmuladd = mat - mmul */
/*inline _Matrix_MUL_ADD_ operator - ( const Matrix& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, -a.alpha, &b, 1, a.t_ab ); }*/

/* mmuladd = mmul + mat_scaled */
inline _Matrix_MUL_ADD_ operator + ( const _Matrix_MUL_& a, const _Matrix_SCALE_& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, b.a, b.alpha, a.t_ab ); }

/* mmuladd = mat_scaled + mmul */
inline _Matrix_MUL_ADD_ operator + ( const _Matrix_SCALE_& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, b.a, b.alpha, a.t_ab ); }

/* mmuladd = mmul - mat_scaled */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_MUL_& a, const _Matrix_SCALE_& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, b.a, -b.alpha, a.t_ab ); }

/* mmuladd = mat_scaled - mmul */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_SCALE_& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, -a.alpha, b.a, b.alpha, a.t_ab ); }

/* mmuladd = mmul + (mat^t) */
inline _Matrix_MUL_ADD_ operator + ( const _Matrix_MUL_& a, const _Matrix_T_& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b.a, b.alpha, a.t_ab + 4 );  }

/* mmuladd = (mat^t) + mmul */
inline _Matrix_MUL_ADD_ operator + ( const _Matrix_T_& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b.a, b.alpha, a.t_ab + 4 );  }

/* mmuladd = mmul - (mat^t) */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_MUL_& a, const _Matrix_T_& b )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha, &b.a, -b.alpha, a.t_ab + 4 );  }

/* mmuladd = (mat^t) - mmul */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_T_& b, const _Matrix_MUL_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, -a.alpha, &b.a, b.alpha, a.t_ab + 4 );  }


/* mmuladd = mat_scaled_shited * mat */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_SCALE_SHIFT_& a, const Matrix& b )
{   return _Matrix_MUL_ADD_( a.a, &b, a.alpha, &b, a.beta, 0 );  }

/* mmuladd = mat * mat_scaled_shited */
inline _Matrix_MUL_ADD_ operator * ( const Matrix& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_MUL_ADD_( &b, a.a, a.alpha, &b, a.beta, 0 );  }

/* mmuladd = mat_scaled_shited * mat_scaled */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_SCALE_SHIFT_& a, const _Matrix_SCALE_& b )
{   return _Matrix_MUL_ADD_( a.a, b.a, a.alpha*b.alpha, b.a, a.beta*b.alpha, 0 );  }

/* mmuladd = mat_scaled * mat_scaled_shited */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_SCALE_& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_MUL_ADD_( b.a, a.a, a.alpha*b.alpha, b.a, a.beta*b.alpha, 0 );  }

/* mmuladd = mat_scaled_shited * (mat^t) */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_SCALE_SHIFT_& a, const _Matrix_T_& b )
{   return _Matrix_MUL_ADD_( a.a, &b.a, a.alpha*b.alpha, &b.a, a.beta*b.alpha, 6 );  }

/* mmuladd = (mat^t) * mat_scaled_shited */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_T_& b, const _Matrix_SCALE_SHIFT_& a )
{   return _Matrix_MUL_ADD_( &b.a, a.a, a.alpha*b.alpha, &b.a, a.beta*b.alpha, 5 );  }

/* mmuladd2 = mmuladd1 * scalar */
inline _Matrix_MUL_ADD_ operator * ( const _Matrix_MUL_ADD_& a, double alpha )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha*alpha, a.c, a.beta*alpha, a.t_abc ); }

/* mmuladd2 = scalar * mmuladd1 */
inline _Matrix_MUL_ADD_ operator * ( double alpha, const _Matrix_MUL_ADD_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, a.alpha*alpha, a.c, a.beta*alpha, a.t_abc ); }

/* -mmuladd */
inline _Matrix_MUL_ADD_ operator - ( const _Matrix_MUL_ADD_& a )
{   return _Matrix_MUL_ADD_( a.a, a.b, -a.alpha, a.c, -a.beta, a.t_abc ); }

/* inv(a)*b, i.e. solve a*x = b */
inline _Matrix_SOLVE_ operator * ( const _Matrix_INV_& a, const Matrix& b )
{   return _Matrix_SOLVE_( &a.a, &b, a.method );    }


/*
* PART III. Logical operations
*/
inline _Matrix_NOT_ operator ~ ( const Matrix& a )
{   return _Matrix_NOT_(&a); }

inline _Matrix_LOGIC_ operator & ( const Matrix& a, const Matrix& b )
{   return _Matrix_LOGIC_( &a, &b, _Matrix_LOGIC_::AND, 0 ); }

inline _Matrix_LOGIC_ operator & ( const _Matrix_NOT_& a, const Matrix& b )
{   return _Matrix_LOGIC_( a.a, &b, _Matrix_LOGIC_::AND, 1 ); }

inline _Matrix_LOGIC_ operator & ( const Matrix& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( &a, b.a, _Matrix_LOGIC_::AND, 2 ); }

inline _Matrix_LOGIC_ operator & ( const _Matrix_NOT_& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( a.a, b.a, _Matrix_LOGIC_::AND, 3 ); }


inline _Matrix_LOGIC_ operator | ( const Matrix& a, const Matrix& b )
{   return _Matrix_LOGIC_( &a, &b, _Matrix_LOGIC_::OR, 0 ); }

inline _Matrix_LOGIC_ operator | ( const _Matrix_NOT_& a, const Matrix& b )
{   return _Matrix_LOGIC_( a.a, &b, _Matrix_LOGIC_::OR, 1 ); }

inline _Matrix_LOGIC_ operator | ( const Matrix& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( &a, b.a, _Matrix_LOGIC_::OR, 2 ); }

inline _Matrix_LOGIC_ operator | ( const _Matrix_NOT_& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( a.a, b.a, _Matrix_LOGIC_::OR, 3 ); }


inline _Matrix_LOGIC_ operator ^ ( const Matrix& a, const Matrix& b )
{   return _Matrix_LOGIC_( &a, &b, _Matrix_LOGIC_::XOR, 0 ); }

inline _Matrix_LOGIC_ operator ^ ( const _Matrix_NOT_& a, const Matrix& b )
{   return _Matrix_LOGIC_( a.a, &b, _Matrix_LOGIC_::XOR, 1 ); }

inline _Matrix_LOGIC_ operator ^ ( const Matrix& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( &a, b.a, _Matrix_LOGIC_::XOR, 2 ); }

inline _Matrix_LOGIC_ operator ^ ( const _Matrix_NOT_& a, const _Matrix_NOT_& b )
{   return _Matrix_LOGIC_( a.a, b.a, _Matrix_LOGIC_::XOR, 3 ); }


inline _Matrix_UN_LOGIC_ operator & ( const Matrix& a, double alpha )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::AND, 0 ); }

inline _Matrix_UN_LOGIC_ operator & ( double alpha, const Matrix& a )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::AND, 0 ); }

inline _Matrix_UN_LOGIC_ operator & ( const _Matrix_NOT_& a, double alpha )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::AND, 1 ); }

inline _Matrix_UN_LOGIC_ operator & ( double alpha, const _Matrix_NOT_& a )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::AND, 1 ); }


inline _Matrix_UN_LOGIC_ operator | ( const Matrix& a, double alpha )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::OR, 0 ); }

inline _Matrix_UN_LOGIC_ operator | ( double alpha, const Matrix& a )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::OR, 0 ); }

inline _Matrix_UN_LOGIC_ operator | ( const _Matrix_NOT_& a, double alpha )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::OR, 1 ); }

inline _Matrix_UN_LOGIC_ operator | ( double alpha, const _Matrix_NOT_& a )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::OR, 1 ); }


inline _Matrix_UN_LOGIC_ operator ^ ( const Matrix& a, double alpha )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::XOR, 0 ); }

inline _Matrix_UN_LOGIC_ operator ^ ( double alpha, const Matrix& a )
{   return _Matrix_UN_LOGIC_( &a, alpha, _Matrix_LOGIC_::XOR, 0 ); }

inline _Matrix_UN_LOGIC_ operator ^ ( const _Matrix_NOT_& a, double alpha )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::XOR, 1 ); }

inline _Matrix_UN_LOGIC_ operator ^ ( double alpha, const _Matrix_NOT_& a )
{   return _Matrix_UN_LOGIC_( a.a, alpha, _Matrix_LOGIC_::XOR, 1 ); }


/*
* PART IV. Comparison operations
*/
inline _Matrix_CMP_ operator > ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_GT ); }

inline _Matrix_CMP_ operator >= ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_GE ); }

inline _Matrix_CMP_ operator < ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_LT ); }

inline _Matrix_CMP_ operator <= ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_LE ); }

inline _Matrix_CMP_ operator == ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_EQ ); }

inline _Matrix_CMP_ operator != ( const Matrix& a, const Matrix& b )
{   return _Matrix_CMP_( &a, &b, CV_CMP_NE ); }


inline _Matrix_CMP_ operator > ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_GT ); }

inline _Matrix_CMP_ operator > ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_LT ); }

inline _Matrix_CMP_ operator >= ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_GE ); }

inline _Matrix_CMP_ operator >= ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_LE ); }

inline _Matrix_CMP_ operator < ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_LT ); }

inline _Matrix_CMP_ operator < ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_GT ); }

inline _Matrix_CMP_ operator <= ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_LE ); }

inline _Matrix_CMP_ operator <= ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_GE ); }

inline _Matrix_CMP_ operator == ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_EQ ); }

inline _Matrix_CMP_ operator == ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_EQ ); }

inline _Matrix_CMP_ operator != ( const Matrix& a, double alpha )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_NE ); }

inline _Matrix_CMP_ operator != ( double alpha, const Matrix& a )
{   return _Matrix_CMP_( &a, alpha, CV_CMP_NE ); }


/*
* PART V. Speedup for some augmented assignments to Matrix
*/

inline Matrix& Matrix::operator += ( const _Matrix_SCALE_& scale_mat )
{   return (*this = *this + scale_mat); }

inline Matrix& Matrix::operator += ( const _Matrix_SCALE_SHIFT_& scale_mat )
{   return (*this = *this + scale_mat); }

inline Matrix& Matrix::operator += ( const _Matrix_MUL_& mmul )
{   return (*this = mmul + *this);  }

inline Matrix& Matrix::operator -= ( const _Matrix_SCALE_& scale_mat )
{   return (*this = *this - scale_mat);  }

inline Matrix& Matrix::operator -= ( const _Matrix_SCALE_SHIFT_& scale_mat )
{   return (*this = *this - scale_mat);  }

inline Matrix& Matrix::operator -= ( const _Matrix_MUL_& mmul )
{   return (*this = -mmul + *this);  }

inline Matrix& Matrix::operator *= ( const _Matrix_SCALE_& scale_mat )
{   return (*this = *this * scale_mat);  }

inline Matrix& Matrix::operator *= ( const _Matrix_SCALE_SHIFT_& scale_mat )
{   return (*this = *this * scale_mat);  }

/****************************************************************************************\
*                        misc. operations on temporary matrices (+,-,*)                  *
\****************************************************************************************/

inline _Matrix_TMP_::_Matrix_TMP_() : Matrix() {}

inline _Matrix_TMP_::_Matrix_TMP_( const _Matrix_TMP_& mat ) : Matrix( (const CvMat&)mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline _Matrix_TMP_& _Matrix_TMP_::operator = ( const _Matrix_TMP_& mat )
{
    if( this != &mat )
    {
        release();
        memcpy( this, &mat, sizeof(mat));
        if( refcount )
            (*refcount)++;
    }
    return *this;
}

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_T_& mat_t ) : Matrix( mat_t )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_INV_& inv_mat ) : Matrix( inv_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_ADD_& mat_add ) : Matrix( mat_add )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_ADD_EX_& mat_add ) : Matrix( mat_add )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_SCALE_& scale_mat ) : Matrix( scale_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_SCALE_SHIFT_& scale_shift_mat ) : Matrix( scale_shift_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_MUL_& mmul ) : Matrix( mmul )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_MUL_ADD_& mmuladd ) : Matrix( mmuladd )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_LOGIC_& mat_logic ) : Matrix( mat_logic )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_UN_LOGIC_& mat_logic ) : Matrix( mat_logic )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_NOT_& not_mat ) : Matrix( not_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_DOT_OP_& dot_mul ) : Matrix( dot_mul )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_SOLVE_& solve_mat ) : Matrix( solve_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_CMP_& cmp_mat ) : Matrix( cmp_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_TMP_& mat, CvRect rect ) :
    Matrix( (const Matrix&)mat, rect )
{   type |= CV_MAT_TEMP_FLAG;   }

inline  _Matrix_TMP_::_Matrix_TMP_( const _Matrix_CVT_& cvt_mat ) : Matrix( cvt_mat )
{   type |= CV_MAT_TEMP_FLAG;   }

inline _Matrix_TMP_::_Matrix_TMP_( const _Matrix_TMP_& mat, int k, int i ) :
    Matrix( (const Matrix&)mat, k, i )
{   type |= CV_MAT_TEMP_FLAG;   }


/*
* the base proxy class implementation
*/

/* a.*b */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::mul( const Matrix& a ) const
{   return ((_Matrix_TMP_)*this).mul(a);  }

/* a.*b*alpha */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::mul( const _Matrix_SCALE_& a ) const
{   return ((_Matrix_TMP_)*this).mul(a);  }

/* a./b */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::div( const Matrix& a ) const
{   return ((_Matrix_TMP_)*this).div(a);  }

/* a./(b*alpha) */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::div( const _Matrix_SCALE_& a ) const
{   return ((_Matrix_TMP_)*this).div(a);  }

/* a.max(b) */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::min( const Matrix& a ) const
{   return ((_Matrix_TMP_)*this).min(a);  }

/* a.min(b) */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::max( const Matrix& a ) const
{   return ((_Matrix_TMP_)*this).max(a);  }

/* a.max(alpha) */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::min( double alpha ) const
{   return ((_Matrix_TMP_)*this).min(alpha);  }

/* a.min(alpha) */
inline _Matrix_DOT_OP_ _Matrix_BASE_OP_::max( double alpha ) const
{   return ((_Matrix_TMP_)*this).max(alpha);  }


inline _Matrix_INV_  _Matrix_BASE_OP_::inv( int method ) const
{   return ((_Matrix_TMP_)*this).inv(method);  }

inline _Matrix_T_  _Matrix_BASE_OP_::t() const
{   return ((_Matrix_TMP_)*this).t();  }

inline _Matrix_CVT_ _Matrix_BASE_OP_::cvt( int newdepth, double scale, double shift ) const
{   return ((_Matrix_TMP_)*this).cvt( newdepth, scale, shift ); }

inline _Matrix_TMP_  _Matrix_BASE_OP_::row( int r ) const
{   return _Matrix_TMP_((_Matrix_TMP_)*this, 0, r ); }

inline _Matrix_TMP_  _Matrix_BASE_OP_::rowrange( int row1, int row2 ) const
{   
    _Matrix_TMP_ m = (_Matrix_TMP_)*this;
    assert( 0 <= row1 && row1 < row2 && row2 <= m.height );
    return _Matrix_TMP_( m, cvRect( 0, row1, m.width, row2 - row1 ));
}

inline _Matrix_TMP_  _Matrix_BASE_OP_::col( int c ) const
{   return _Matrix_TMP_( (_Matrix_TMP_)*this, 1, c ); }

inline _Matrix_TMP_  _Matrix_BASE_OP_::colrange( int col1, int col2 ) const
{   
    _Matrix_TMP_ m = (_Matrix_TMP_)*this;
    assert( 0 <= col1 && col1 < col2 && col2 <= m.width );
    return _Matrix_TMP_( m, cvRect( col1, 0, col2 - col1, m.height ));
}

inline _Matrix_TMP_  _Matrix_BASE_OP_::rect( CvRect r ) const
{   return _Matrix_TMP_( (_Matrix_TMP_)*this, r ); }

inline _Matrix_TMP_  _Matrix_BASE_OP_::diag( int d ) const
{   return _Matrix_TMP_( (_Matrix_TMP_)*this, -1, d ); }

inline double _Matrix_BASE_OP_::det() const
{   return ((_Matrix_TMP_)*this).det(); }

inline double _Matrix_BASE_OP_::norm( int norm_type ) const
{   return ((_Matrix_TMP_)*this).norm( norm_type ); }

inline CvScalar _Matrix_BASE_OP_::sum() const
{   return ((_Matrix_TMP_)*this).sum(); }

inline double _Matrix_BASE_OP_::min( CvPoint* minloc ) const
{   return ((_Matrix_TMP_)*this).min( minloc ); }

inline double _Matrix_BASE_OP_::max( CvPoint* maxloc ) const
{   return ((_Matrix_TMP_)*this).max( maxloc ); }


/****************************************************************************************/
/*                              proxy classes implementation.                           */
/*                              part I. constructors                                    */
/****************************************************************************************/

/* constructors */
inline _Matrix_COPY_::_Matrix_COPY_( const Matrix* _a ) : a((Matrix*)_a)  {}

inline _Matrix_CVT_::_Matrix_CVT_( const Matrix* _a, int _newdepth,
                                 double _scale, double _shift ) :
    a(*(Matrix*)_a), newdepth(_newdepth), scale(_scale), shift(_shift)  {}

inline _Matrix_T_::_Matrix_T_( const Matrix* _a ) : a(*(Matrix*)_a), alpha(1)  {}


inline _Matrix_T_::_Matrix_T_( const Matrix* _a, double _alpha ) :
    a(*(Matrix*)_a), alpha(_alpha)  {}


inline _Matrix_INV_::_Matrix_INV_( const Matrix* _a, int _method ) :
    a(*(Matrix*)_a), method(_method) {}


inline _Matrix_MUL_::_Matrix_MUL_( const Matrix* _a, const Matrix* _b, int _t_ab ) :
    a((Matrix*)_a), b((Matrix*)_b), alpha(1), t_ab(_t_ab) {}


inline _Matrix_MUL_::_Matrix_MUL_( const Matrix* _a, const Matrix* _b,
                                 double _alpha, int _t_ab ) :
    a((Matrix*)_a), b((Matrix*)_b), alpha(_alpha), t_ab(_t_ab) {}


inline _Matrix_MUL_ADD_::_Matrix_MUL_ADD_( const Matrix* _a, const Matrix* _b,
                                         const Matrix* _c, int _t_abc ) :
    a((Matrix*)_a), b((Matrix*)_b), c((Matrix*)_c), t_abc(_t_abc) {}


inline _Matrix_MUL_ADD_::_Matrix_MUL_ADD_( const Matrix* _a, const Matrix* _b, double _alpha,
                                         const Matrix* _c, double _beta, int _t_abc ) :
    a((Matrix*)_a), b((Matrix*)_b), alpha(_alpha),
    c((Matrix*)_c), beta(_beta), t_abc(_t_abc) {}


inline _Matrix_ADD_::_Matrix_ADD_( const Matrix* _a, const Matrix* _b, double _beta ) :
    a((Matrix*)_a), b((Matrix*)_b), beta(_beta) {}


inline _Matrix_ADD_EX_::_Matrix_ADD_EX_( const Matrix* _a, double _alpha,
                                       const Matrix* _b, double _beta, double _gamma ) :
    a((Matrix*)_a), alpha(_alpha), b((Matrix*)_b), beta(_beta), gamma(_gamma) {}


inline _Matrix_SCALE_::_Matrix_SCALE_( const Matrix* _a, double _alpha ) :
    a((Matrix*)_a), alpha(_alpha) {}


inline _Matrix_SCALE_SHIFT_::_Matrix_SCALE_SHIFT_( const Matrix* _a,
                                                 double _alpha, double _beta ) :
    a((Matrix*)_a), alpha(_alpha), beta(_beta) {}


inline _Matrix_LOGIC_::_Matrix_LOGIC_( const Matrix* _a, const Matrix* _b,
                                            _Matrix_LOGIC_::Op _op, int _flags ) :
    a((Matrix*)_a), b((Matrix*)_b), op(_op), flags(_flags) {}


inline _Matrix_UN_LOGIC_::_Matrix_UN_LOGIC_( const Matrix* _a, double _alpha,
                                           _Matrix_LOGIC_::Op _op, int _flags ) :
    a((Matrix*)_a), alpha(_alpha), op(_op), flags(_flags) {}


inline _Matrix_NOT_::_Matrix_NOT_( const Matrix* _a ) :
    a((Matrix*)_a) {}


inline _Matrix_DOT_OP_::_Matrix_DOT_OP_( const Matrix* _a, const Matrix* _b,
                                       int _op, double _alpha ) :
    a(*_a), b((Matrix*)_b), op(_op), alpha(_alpha) {}


inline _Matrix_SOLVE_::_Matrix_SOLVE_( const Matrix* _a, const Matrix* _b, int _method ) :
    a((Matrix*)_a), b((Matrix*)_b), method(_method) {}

inline _Matrix_CMP_::_Matrix_CMP_( const Matrix* _a, const Matrix* _b, int _cmp_op ) :
    a((Matrix*)_a), b((Matrix*)_b), alpha(0), cmp_op(_cmp_op) {}

inline _Matrix_CMP_::_Matrix_CMP_( const Matrix* _a, double _alpha, int _cmp_op ) :
    a((Matrix*)_a), b(0), alpha(_alpha), cmp_op(_cmp_op) {}

/****************************************************************************************/
/*                              proxy classes implementation.                           */
/*                              part II. conversion to _Matrix_TMP_                      */
/****************************************************************************************/

inline _Matrix_T_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_INV_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_MUL_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_SCALE_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_SCALE_SHIFT_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_ADD_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_ADD_EX_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_MUL_ADD_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_LOGIC_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_UN_LOGIC_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_NOT_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_DOT_OP_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_SOLVE_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_CMP_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_( *this );    }

inline _Matrix_CVT_::operator _Matrix_TMP_() const
{   return _Matrix_TMP_(*this);   }

inline _Matrix_COPY_::operator Matrix() const
{   return *a;   }

/****************************************************************************************/
/*                              proxy classes implementation.                           */
/*                              part III. custom overrided methods                      */
/****************************************************************************************/

inline _Matrix_DOT_OP_ _Matrix_SCALE_::mul( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( a, &mat, '*', alpha );   }

inline _Matrix_DOT_OP_ _Matrix_SCALE_::mul( const _Matrix_SCALE_& mat ) const
{   return _Matrix_DOT_OP_( a, mat.a, '*', alpha*mat.alpha );   }

inline _Matrix_DOT_OP_ _Matrix_SCALE_::div( const Matrix& mat ) const
{   return _Matrix_DOT_OP_( a, &mat, '/', alpha );   }

inline _Matrix_DOT_OP_ _Matrix_SCALE_::div( const _Matrix_SCALE_& mat ) const
{   return _Matrix_DOT_OP_( a, mat.a, '/', alpha/mat.alpha );   }

inline _Matrix_DOT_OP_ operator * ( const _Matrix_DOT_OP_& dot_op, double alpha )
{   return _Matrix_DOT_OP_( &dot_op.a, dot_op.b, dot_op.op, dot_op.alpha * alpha );  }

inline _Matrix_DOT_OP_ operator * ( double alpha, const _Matrix_DOT_OP_& dot_op )
{   return _Matrix_DOT_OP_( &dot_op.a, dot_op.b, dot_op.op, dot_op.alpha * alpha );  }

inline _Matrix_DOT_OP_ operator / ( double alpha, const Matrix& mat )
{   return _Matrix_DOT_OP_( &mat, 0, '/', alpha );  }

inline _Matrix_DOT_OP_ operator / ( double alpha, const _Matrix_SCALE_& mat )
{   return _Matrix_DOT_OP_( mat.a, 0, '/', alpha/mat.alpha );  }


inline double _Matrix_T_::det() const
{   return a.det();     }

inline double _Matrix_T_::norm( int norm_type ) const
{   return a.norm( norm_type );    }

inline double _Matrix_ADD_::norm( int norm_type ) const
{
    if( beta == -1 )
        return cvNorm( a, b, norm_type );
    else
        return ((_Matrix_TMP_)*this).norm( norm_type );
}

inline _Matrix_DOT_OP_ _Matrix_ADD_::abs() const
{
    if( beta == -1 )
        return _Matrix_DOT_OP_( a, b, 'a', 0 );
    else
        return ((_Matrix_TMP_)*this).abs();
}

inline _Matrix_DOT_OP_ _Matrix_SCALE_SHIFT_::abs() const
{
    if( alpha == 1 )
        return _Matrix_DOT_OP_( a, 0, 'a', -beta );
    else
        return ((_Matrix_TMP_)*this).abs();
}

inline std::ostream & print_row(std::ostream & out, const Matrix &m, const int & i){
	out<<"[";
	for(int j=0; j<(m.cols-1); j++){
		out<<cvmGet(&m,i,j)<<", ";
	}
	out<<cvmGet(&m,i,m.cols-1)<<"]";
	return out;
}
inline std::ostream & operator <<(std::ostream & out, const Matrix &m) {
	out<<"[";
	for(int i=0; i<(m.rows-1); i++){
		print_row(out, m, i);
		out<<",\n";
	}
	print_row(out, m, m.rows-1);
	out<<"]";
	return out;
}
#endif /*_CVMAT_HPP_*/

