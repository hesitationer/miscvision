#ifndef __GAUSSIAN_ND_HH
#define __GAUSSIAN_ND_HH

#include <stdlib.h> //drand48
#include <math.h> //M_SQRT1_2 , M_PI, tan2
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <cv/matrix.hh>
#include <cv/gaussian.h>

class GaussianND { 
protected:
	Matrix<float> m_mu;
    Matrix<float> m_covar;
    Matrix<float> m_inverse;
    Matrix<float> m_cholesky;
    double m_det;
public:	
	GaussianND(int size): m_mu(size,1), m_covar(size,size), m_inverse(size, size), m_cholesky(size, size){
    }
    
    GaussianND(const GaussianND & g): m_mu(g.m_mu), m_covar(g.m_covar), 
									  m_inverse(g.m_inverse), m_cholesky(g.m_cholesky), m_det(g.m_det)
	{
	}	
	GaussianND(const char * fname ){
		this->load(fname);
	}
    
	// pdf = exp( -.5 * x^T*inv(K)*x) / sqrt( (2*PI)^d * det(K) )
    float pdf(const Matrix<float> & x){
		Matrix<float> z(m_mu.rows, m_mu.cols);
		Matrix<float> zi(m_mu.rows, m_mu.cols);
		Matrix<float> evals(m_mu.rows, m_mu.cols);
		double xtx;


		z = x - m_mu;
		zi = m_inverse * z;
		xtx = cvDotProduct( &zi, &z );
		cvCheckArr(&z);
		cvCheckArr(&zi);

        //cvSVD( &m_covar, &evals);
		//std::cout<<evals<<std::endl;
		//std::cout<<"xtx="<<xtx<<std::endl;
		//std::cout<<"det="<<m_det<<std::endl;
		double norm = fabs(m_det);
		norm *= pow(2*M_PI, m_mu.rows);
		//std::cout<<"norm="<<sqrt(norm)<<std::endl;
		double p = exp(-.5 * xtx) / sqrt(norm);
		//std::cout<<"gauss="<<p<<std::endl;
		//std::cout<<"exp="<<exp(-.5*xtx)<<std::endl;
        
        return p;
    }
	
	static double x_Sigma_x(const CvArr * x, const CvArr * mean, const CvArr * covar_inv, 
			CvArr * tmp1, CvArr * tmp2){
		CvMat tmp1_stub;
		CvMat tmp2_stub;
		cvSub(x, mean, tmp1);

		tmp1 = cvReshape( tmp1, &tmp1_stub, 1, 1);
		tmp2 = cvReshape( tmp2, &tmp2_stub, 1, 1);
		cvMatMul( tmp1, covar_inv, tmp2 );
		return cvDotProduct( tmp1, tmp2 );
	}

	static double log_likelihood_ratio(const Matrix<float> &x, const GaussianND & g1, const GaussianND & g2 ){
		Matrix<float> tmp1(g1.mu().rows, g1.mu().cols);
		Matrix<float> tmp2(g1.mu().rows, g1.mu().cols);
		double offset = log(sqrt( fabs( g1.det() ) )) - log(sqrt( fabs( g2.det() ) ));
		double llr = GaussianND::x_Sigma_x( &x, &g1.mu(), &g1.sigma_inv(), &tmp1, &tmp2) -
			   	     GaussianND::x_Sigma_x( &x, &g2.mu(), &g2.sigma_inv(), &tmp1, &tmp2) - offset;
		return llr;
	}

	// decision function
	float g(const Matrix<float> & x){
		Matrix<float> z(m_mu.rows, m_mu.cols);
		Matrix<float> zi(m_mu.rows, m_mu.cols);
		Matrix<float> evals(m_mu.rows, m_mu.cols);
		double xtx;


		z = x - m_mu;
		zi = m_inverse * z;
		xtx = cvDotProduct( &zi, &z );

		return -.5 * (xtx + log( m_det ));	
	}

    Matrix<float> sample() const{
        Matrix<float> x(m_mu.rows, m_mu.cols);
        Gaussian n(0,1);
        for(int i=0; i<x.rows; i++){
            x[i][0] =  n.sample();
        }
       	cvMul(&m_cholesky, &x, &x);
		cvAdd(&x, &m_mu, &x);
		return x;
    }
	void setMean(const Matrix<float> &x){
        cvCopy(&x, &m_mu);
	}
	const Matrix<float> & mu() const{
		return m_mu;
	}
	Matrix<float> & mu() {
		return m_mu;
	}
	const Matrix<float> & sigma() const{
		return m_covar;
	}
	Matrix<float> & sigma(){
		return m_covar;
	}
	const Matrix<float> & sigma_inv() const{
		return m_inverse;
	}
	Matrix<float> & sigma_inv(){
		return m_inverse;
	}
	double det() const { 
		return m_det;
	}
	const Matrix<float> &getMean() const{
		return m_mu;
	}
    const Matrix<float> & getCovar() const{
        return m_covar;
    }
    void setCovar(const Matrix<float> & v){
        cvCopy(&v, &m_covar);
		this->recalc();
    }
	void recalc(){
		// for pdf
		cvInvert( &m_covar, &m_inverse, CV_SVD_SYM );
		m_det=cvDet( &m_covar );

        // for sampling
		cvCholesky(&m_covar, &m_cholesky);
	}
	bool save( const char * fname ){
		bool ret = true;
		char text[256];
		if(mkdir( fname, 0777 )!=0 && errno!=EEXIST){
			perror("mkdir");
			return false;
		}
		snprintf(text, 256, "%s/mean.cvm", fname);
		ret = ret && m_mu.save( text );
		snprintf(text, 256, "%s/covar.cvm", fname);
		ret = ret && m_covar.save( text );
		snprintf(text, 256, "%s/covar_inv.cvm", fname);
		ret = ret && m_inverse.save( text );
		snprintf(text, 256, "%s/cholesky.cvm", fname);
		ret = ret && m_inverse.save( text );

		return ret;
	}
	bool load( const char * fname ){
		bool ret = true;
		char text[256];
        snprintf(text, 256, "%s/mean.cvm", fname);
        ret = ret && m_mu.load( text );
        snprintf(text, 256, "%s/covar.cvm", fname);
        ret = ret && m_covar.load( text );
        snprintf(text, 256, "%s/covar_inv.cvm", fname);
        ret = ret && m_inverse.load( text );
        snprintf(text, 256, "%s/cholesky.cvm", fname);
        ret = ret && m_inverse.load( text );
		return ret;
	}
};


#endif //__GAUSSIAN_ND_HH
