#ifndef __GAUSSIAN_HH
#define __GAUSSIAN_HH

#include <stdlib.h> //drand48
#include <math.h> //M_SQRT1_2 , M_PI, tan2
#include <assert.h>
#include <stdio.h>
#include <ostream>
#include <cv/wincompat.h> //drand48, etc for windows
using namespace std;

#define DEG2RAD(x) (x*M_PI/180.0)
#define RAD2DEG(x) (x*180*M_1_PI)

//generate a gaussian sample - 
//courtesy of http://www-math.mit.edu/~spielman/ECC/randGen.c
static float gasdev()
{
        int iset=0;
        float gset;
        float fac,rsq,v1,v2;

        if (iset == 0) {
                do {
                        v1=2.0*drand48()-1.0;
                        v2=2.0*drand48()-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq==0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } 
        else {
                iset =0;
                return gset;
        }
        return 0;
}

class Gaussian {
protected:
	float m_mu;
	float m_sig;

public:	
	Gaussian(float m=0, float s=1): m_mu(m), m_sig(s){}
    Gaussian(const Gaussian & g): m_mu(g.m_mu), m_sig(g.m_sig){}	
	void pdf(IplImage * src, IplImage * res){
		this->z(src, res);
		cvMul(res, res, res, -0.5);
		cvExp(res, res);
		cvScale(res, res, sqrt(2*M_PI*m_sig*m_sig));
	}
    double pdf(double x){
        return exp(-.5 * z(x)*z(x)) / sqrt(2*M_PI*m_sig*m_sig);
    }
    double phi(double z){
		return .5 * erf(z * M_SQRT1_2);
	}
	void z(IplImage * src, IplImage * res){
		cvSubS(src, cvScalar(m_mu), res);
		cvScale(res, res, 1/(m_sig*m_sig));
	}
	double z(double x){
		return (x-m_mu) / (m_sig*m_sig);
	}
	double cdf(double l, double u){
		return phi(z(u)) - phi(z(l));
	}
	void mul(const Gaussian & g){
		m_mu = m_mu/this->getVar() + g.m_mu/g.getVar();
		m_sig = 1.0/this->getVar() + 1.0/g.getVar();
		m_mu /= m_sig;
		m_sig = 1.0/sqrt(m_sig);
	}
    double sample(){
        return m_mu + gasdev()*m_sig*m_sig;
    }
	void setMean(float x){
		m_mu = x;
	}
	float getMean() const{
		return m_mu;
	}
	float getVar() const{
		return m_sig*m_sig;
	}
	float getSig() const{
		return m_sig;
	}
	void setSig(float x){
		m_sig = x;
	}
	void print(){
		printf("N(%f, %f)\n", m_mu, m_sig);
	}
    friend ostream & operator<< (ostream &s, const Gaussian &g){
        return s<<"N("<<g.m_mu<<", "<<g.m_sig<<")";
    }
};


#endif //__STAT_HH
