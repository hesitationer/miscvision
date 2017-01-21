
// this file is included by matrix.hpp and is not meant to be included seperately

namespace cv {
	
/// base Matrix operation class
template <typename pixel_t,int nch> struct Matrix_OP {
	virtual void assign(Matrix<pixel_t, nch> & X) const = 0;
};

/// Generalized Matrix Multiplucation operations
template <typename pixel_t, int nch> struct Matrix_GEMM : Matrix_OP<pixel_t, nch> {
	typedef Matrix<pixel_t,nch> array_t;
	const array_t * m_A;
	const array_t * m_B;
	const array_t * m_C;
	double _alpha;
	double _beta;
	int _tABC;
	Matrix_GEMM(const array_t *A, const array_t *B, double alpha, const array_t * C, double beta, int tABC=0):
		m_A(A),
		m_B(B),
		m_C(C),
		_alpha(alpha),
		_beta(beta),
		_tABC(tABC)
	{
	}
	virtual void assign(Matrix<pixel_t,nch> & X) const{
		if(m_B!=NULL){
			cvGEMM(m_A, m_B, _alpha, m_C, _beta, &X, _tABC);
		}
	}
};

/// Add scalar to matrix
template <typename pixel_t, int nch> struct Matrix_ADDS : Matrix_OP<pixel_t, nch>{
	const Matrix<pixel_t,nch> * m_A;
	CvScalar _shift;
	Matrix_ADDS(const Matrix<pixel_t,nch> *A, CvScalar shift=0):
		        m_A(A),
				_shift(shift)
	{
	}
	virtual void assign(Matrix<pixel_t,nch> & X) const{
		cvAddS(m_A, _shift, &X);
	}

};

/// Add two matrices
template <typename pixel_t, int nch> struct Matrix_ADD : Matrix_OP<pixel_t, nch>{
	const Matrix<pixel_t,nch> * m_A;
	const Matrix<pixel_t,nch> * m_B;
	double _scale1;
	double _scale2;
	double _shift;
	Matrix_ADD(const Matrix<pixel_t,nch> *A, const Matrix<pixel_t,nch> *B, double scale1, double scale2=1, double shift=0):
		m_A(A),
		m_B(B),
		_scale1(scale1),
		_scale2(scale2),
		_shift(shift)

	{
	}
	virtual void assign(Matrix<pixel_t,nch> & X) const{
		cvAddWeighted(m_A, _scale1, m_B, _scale2, _shift, &X);
	}
};

/// Multiply matrix by a scalar
template <typename pixel_t, int nch> struct Matrix_SCALE : Matrix_OP<pixel_t, nch> {
	const Matrix<pixel_t,nch> * m_A;
	double _scale;
	Matrix_SCALE(const Matrix<pixel_t, nch> * A, double scale):
		m_A(A),
		_scale(scale){}
	virtual void assign(Matrix<pixel_t,nch> & X) const{
		cvScale(m_A, &X, _scale);
	}
};

/// Multiply two matrices
template <typename pixel_t, int nch> struct Matrix_MUL : Matrix_GEMM<pixel_t, nch> {
	typedef Matrix<pixel_t,nch> array_t;
	Matrix_MUL(const array_t *A, const array_t *B):
		Matrix_GEMM<pixel_t,nch>(A,B,1.0,NULL,1.0)
	{
	}
};

// scale * A
template <typename pixel_t, int nch>
Matrix_SCALE<pixel_t,nch> operator * (const double & scale, const Matrix<pixel_t,nch> & A) {
	return Matrix_SCALE<pixel_t,nch>(&A, scale);
}
// A * scale
template <typename pixel_t, int nch>
Matrix_SCALE<pixel_t,nch> operator * (const Matrix<pixel_t,nch> & A, const double & scale) {
	return Matrix_SCALE<pixel_t,nch>(&A, scale);
}

// scale * (A* scale_n) == A*(scale*scale_n)
template <typename pixel_t, int nch>
Matrix_SCALE<pixel_t,nch> operator * (const double & scale, const Matrix_SCALE<pixel_t,nch> & A) {
	return Matrix_SCALE<pixel_t,nch>(A.m_A, A._scale*scale);
}

// (A* scale_n) * scale == A*(scale*scale_n)
template <typename pixel_t, int nch>
Matrix_SCALE<pixel_t,nch> operator * (const Matrix_SCALE<pixel_t,nch> & A, const double & scale) {
	return Matrix_SCALE<pixel_t,nch>(A.m_A, A._scale*scale);
}

// (A * scale) + B
template <typename pixel_t, int nch>
Matrix_ADD<pixel_t,nch> operator + (const Matrix_SCALE<pixel_t,nch> & A, const Matrix<pixel_t,nch> & B){
	return Matrix_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);
}

// A  + b
template <typename pixel_t, int nch>
Matrix_ADDS<pixel_t,nch> operator + (const Matrix<pixel_t,nch> & A, const CvScalar & shift){
	return Matrix_ADDS<pixel_t,nch>(&A, shift);
}

// A  - b
template <typename pixel_t, int nch>
Matrix_ADDS<pixel_t,nch> operator - (const Matrix<pixel_t,nch> & A, const CvScalar & shift){
	return Matrix_ADDS<pixel_t,nch>(&A, cvScalar(-shift.val[0], -shift.val[1], -shift.val[2], -shift.val[3]));
}

// A  + b
template <typename pixel_t, int nch>
Matrix_ADDS<pixel_t,nch> operator + (const Matrix<pixel_t,nch> & A, const double & shift){
	return Matrix_ADDS<pixel_t,nch>(&A, cvScalarAll(shift));
}

// A  - b
template <typename pixel_t, int nch>
Matrix_ADDS<pixel_t,nch> operator - (const Matrix<pixel_t,nch> & A, const double & shift){
	return Matrix_ADDS<pixel_t,nch>(&A, cvScalarAll(-shift));
}

// B + (A * scale)
template <typename pixel_t, int nch>
Matrix_ADD<pixel_t,nch> operator + (const Matrix<pixel_t,nch> & B, const Matrix_SCALE<pixel_t,nch> & A){
	return Matrix_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);
}
// (A*scale) + (B * scale)
template <typename pixel_t, int nch>
Matrix_ADD<pixel_t,nch> operator + (const Matrix_SCALE<pixel_t,nch> & A, const Matrix_SCALE<pixel_t,nch> & B){
	return Matrix_ADD<pixel_t,nch>(A.m_A, B.m_A, A._scale, B._scale);
}

// (A+B) * scale
template <typename pixel_t, int nch>
Matrix_ADD<pixel_t,nch> operator * (const Matrix_ADD<pixel_t,nch> & A, const double & scale){
	return Matrix_ADD<pixel_t,nch>(A.m_A, A.m_B, A._scale1*scale, A._scale2*scale);
}

// A*B
template <typename pixel_t, int nch>
Matrix_MUL<pixel_t,nch> operator * (const Matrix<pixel_t,nch> & A, const Matrix<pixel_t,nch> &B) {
	return Matrix_MUL<pixel_t,nch>(&A, &B);
}

// A+B
template <typename pixel_t,int nch>
Matrix_ADD<pixel_t, nch> operator + (const Matrix<pixel_t,nch> & A, const Matrix<pixel_t,nch> &B) {
	return Matrix_ADD<pixel_t,nch>(&A, &B, 1.0, 1.0);
}

// A-B
template <typename pixel_t, int nch>
Matrix_ADD<pixel_t, nch> operator - (const Matrix<pixel_t, nch> & A, const Matrix<pixel_t, nch> &B) {
	return Matrix_ADD<pixel_t, nch>(&A, &B, 1.0, -1.0);
}

// A*B+C
template <typename pixel_t, int nch>
Matrix_GEMM<pixel_t, nch> operator + (Matrix_MUL<pixel_t, nch> & M, const Matrix<pixel_t, nch> &C){
	return Matrix_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);
}

//C+A*B
template <typename pixel_t, int nch>
Matrix_GEMM<pixel_t, nch> operator + (const Matrix<pixel_t, nch> &C, Matrix_MUL<pixel_t, nch> & M){
	return Matrix_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);
}

//C - A*B
template <typename pixel_t, int nch>
Matrix_GEMM<pixel_t, nch> operator - (const Matrix<pixel_t, nch> &C, Matrix_MUL<pixel_t, nch> & M){
	return Matrix_GEMM<pixel_t, nch>(M.m_A, M.m_B, -1.0, &C, 1.0);
}

//A*B - C
template <typename pixel_t, int nch>
Matrix_GEMM<pixel_t, nch> operator - (Matrix_MUL<pixel_t, nch> & M, const Matrix<pixel_t, nch> &C){
	return Matrix_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, -1.0);
}

} // namespace cv
