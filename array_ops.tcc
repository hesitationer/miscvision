
// this file is included by matrix.hpp and is not meant to be included seperately

namespace cv {
	
/// base cv::Array operation class
template <typename pixel_t, int nch> struct Array_OP {
	virtual void assign(CvArr * X) const = 0;
};

/// Generalized Array Multiplucation operations
template <typename pixel_t, int nch> struct Array_GEMM : Array_OP<pixel_t, nch> {
	const CvArr * m_A;
	const CvArr * m_B;
	const CvArr * m_C;
	double _alpha;
	double _beta;
	int _tABC;
	Array_GEMM(const CvArr *A, const CvArr *B, double alpha, const CvArr * C, double beta, int tABC=0):
		m_A(A),
		m_B(B),
		m_C(C),
		_alpha(alpha),
		_beta(beta),
		_tABC(tABC)
	{
	}
	virtual void assign(CvArr * X) const{
		if(m_B!=NULL){
			cvGEMM(m_A, m_B, _alpha, m_C, _beta, X, _tABC);
		}
	}
};

/// Add scalar to matrix
template <typename pixel_t, int nch> struct Array_ADDS : Array_OP<pixel_t, nch>{
	const CvArr * m_A;
	CvScalar _shift;
	Array_ADDS(const CvArr *A, CvScalar shift=0):
		        m_A(A),
				_shift(shift)
	{
	}
	virtual void assign(CvArr * X) const{
		cvAddS(m_A, _shift, X);
	}
};

/// Add two matrices
template <typename pixel_t, int nch> struct Array_ADD : Array_OP<pixel_t, nch>{
	const CvArr * m_A;
	const CvArr * m_B;
	double _scale1;
	double _scale2;
	double _shift;
	Array_ADD(const CvArr *A, const CvArr *B, double scale1, double scale2=1, double shift=0):
		m_A(A),
		m_B(B),
		_scale1(scale1),
		_scale2(scale2),
		_shift(shift)

	{
	}
	virtual void assign(CvArr * X) const{
		cvAddWeighted(m_A, _scale1, m_B, _scale2, _shift, X);
	}
};

/// Multiply matrix by a scalar
template <typename pixel_t, int nch> struct Array_SCALE : Array_OP<pixel_t, nch> {
	const CvArr * m_A;
	double _scale;
	Array_SCALE(const CvArr * A, double scale):
		m_A(A),
		_scale(scale){}
	virtual void assign(CvArr * X) const{
		cvScale(m_A, &X, _scale);
	}
};

/// Multiply two matrices
template <typename pixel_t, int nch> struct Array_MUL : Array_GEMM<pixel_t, nch> {
	Array_MUL(const CvArr *A, const CvArr *B):
		Array_GEMM<pixel_t,nch>(A,B,1.0,NULL,1.0)
	{
	}
};

// scale * A
#define ARRAY_OPS_IMPL(array_t) \
template <typename pixel_t, int nch> \
cv::Array_SCALE<pixel_t,nch> operator * (const double & scale, const array_t <pixel_t,nch> & A) { \
	return cv::Array_SCALE<pixel_t,nch>(&A, scale); \
} \
\
\
/* A * scale */ \
template <typename pixel_t, int nch> \
cv::Array_SCALE<pixel_t,nch> operator * (const array_t<pixel_t,nch> & A, const double & scale) { \
	return cv::Array_SCALE<pixel_t,nch>(&A, scale); \
} \
\
/* scale * (A* scale_n) == A*(scale*scale_n) */ \
template <typename pixel_t, int nch> \
cv::Array_SCALE<pixel_t,nch> operator * (const double & scale, const cv::Array_SCALE<pixel_t,nch> & A) { \
	return cv::Array_SCALE<pixel_t,nch>(A.m_A, A._scale*scale); \
} \
\
/* (A* scale_n) * scale == A*(scale*scale_n) */ \
template <typename pixel_t, int nch> \
cv::Array_SCALE<pixel_t,nch> operator * (const cv::Array_SCALE<pixel_t,nch> & A, const double & scale) { \
	return cv::Array_SCALE<pixel_t,nch>(A.m_A, A._scale*scale); \
}\
\
/* (A * scale) + B */\
template <typename pixel_t, int nch>\
cv::Array_ADD<pixel_t,nch> operator + (const cv::Array_SCALE<pixel_t,nch> & A, const array_t<pixel_t,nch> & B){\
	return cv::Array_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);\
}\
\
/* A  + b */\
template <typename pixel_t, int nch>\
cv::Array_ADDS<pixel_t,nch> operator + (const array_t<pixel_t,nch> & A, const CvScalar & shift){\
	return cv::Array_ADDS<pixel_t,nch>(&A, shift);\
}\
\
/* A  - b */\
template <typename pixel_t, int nch>\
cv::Array_ADDS<pixel_t,nch> operator - (const array_t<pixel_t,nch> & A, const CvScalar & shift){\
	return cv::Array_ADDS<pixel_t,nch>(&A, cvScalar(-shift.val[0], -shift.val[1], -shift.val[2], -shift.val[3]));\
}\
\
/* A  + b */\
template <typename pixel_t, int nch>\
cv::Array_ADDS<pixel_t,nch> operator + (const array_t<pixel_t,nch> & A, const double & shift){\
	return cv::Array_ADDS<pixel_t,nch>(&A, cvScalarAll(shift));\
}\
\
/* A  - b */\
template <typename pixel_t, int nch>\
cv::Array_ADDS<pixel_t,nch> operator - (const array_t<pixel_t,nch> & A, const double & shift){\
	return cv::Array_ADDS<pixel_t,nch>(&A, cvScalarAll(-shift));\
}\
\
/* B + (A * scale) */\
template <typename pixel_t, int nch>\
cv::Array_ADD<pixel_t,nch> operator + (const array_t<pixel_t,nch> & B, const cv::Array_SCALE<pixel_t,nch> & A){\
	return cv::Array_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);\
}\
/* (A*scale) + (B * scale) */\
template <typename pixel_t, int nch>\
cv::Array_ADD<pixel_t,nch> operator + (const cv::Array_SCALE<pixel_t,nch> & A, const cv::Array_SCALE<pixel_t,nch> & B){\
	return cv::Array_ADD<pixel_t,nch>(A.m_A, B.m_A, A._scale, B._scale);\
}\
\
/* (A+B) * scale */\
template <typename pixel_t, int nch>\
cv::Array_ADD<pixel_t,nch> operator * (const cv::Array_ADD<pixel_t,nch> & A, const double & scale){\
	return cv::Array_ADD<pixel_t,nch>(A.m_A, A.m_B, A._scale1*scale, A._scale2*scale);\
}\
\
/* A*B */\
template <typename pixel_t, int nch>\
cv::Array_MUL<pixel_t,nch> operator * (const array_t<pixel_t,nch> & A, const array_t<pixel_t,nch> &B) {\
	return cv::Array_MUL<pixel_t,nch>(&A, &B);\
}\
\
/* A+B */\
template <typename pixel_t,int nch>\
cv::Array_ADD<pixel_t, nch> operator + (const array_t<pixel_t,nch> & A, const array_t<pixel_t,nch> &B) {\
	return cv::Array_ADD<pixel_t,nch>(&A, &B, 1.0, 1.0);\
}\
\
/* A-B */\
template <typename pixel_t, int nch>\
cv::Array_ADD<pixel_t, nch> operator - (const array_t<pixel_t, nch> & A, const array_t<pixel_t, nch> &B) {\
	return cv::Array_ADD<pixel_t, nch>(&A, &B, 1.0, -1.0);\
}\
\
/* A*B+C */\
template <typename pixel_t, int nch>\
cv::Array_GEMM<pixel_t, nch> operator + (cv::Array_MUL<pixel_t, nch> & M, const array_t<pixel_t, nch> &C){\
	return cv::Array_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);\
}\
\
/*C+A*B */\
template <typename pixel_t, int nch>\
cv::Array_GEMM<pixel_t, nch> operator + (const array_t<pixel_t, nch> &C, cv::Array_MUL<pixel_t, nch> & M){\
	return cv::Array_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);\
}\
\
/*C - A*B */\
template <typename pixel_t, int nch>\
cv::Array_GEMM<pixel_t, nch> operator - (const array_t<pixel_t, nch> &C, cv::Array_MUL<pixel_t, nch> & M){\
	return cv::Array_GEMM<pixel_t, nch>(M.m_A, M.m_B, -1.0, &C, 1.0);\
}\
\
/*A*B - C */\
template <typename pixel_t, int nch>\
cv::Array_GEMM<pixel_t, nch> operator - (cv::Array_MUL<pixel_t, nch> & M, const array_t<pixel_t, nch> &C){\
	return cv::Array_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, -1.0);\
}

} // namespace cv
