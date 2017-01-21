#ifndef __VECTOR_HH
#define __VECTOR_HH

#include <math.h>
#include <cv.h>
#include <stdio.h>

#include "matrix.hh"

template <class T>
class Vector: public Matrix<T> {
    
public:
	Vector() : Matrix<T>(1,1){}
    Vector(const Vector<T> & v){
		(*this) = v;
	}
    Vector(int length) : Matrix<T>(length, 1){}
    Vector(const Matrix<T> & m){
		(*this) = m;
	}
    const T & get(const int &i) const{
        return (*this)(i);
    }
    void set(const int &i, const T & t){
        (*this)(i) = t;
    }
    const T & operator[] (const int &i) const{
        #ifdef DEBUG
            assert(i<this->rows);
            assert(i>=0);
            assert(cvPtr2D(this,i,0)==this->data.ptr+i*this->step);
        #endif
        return ((T*)(this->data.ptr+i*this->step))[0];
    }
    const T & operator() (const int &i) const{
        #ifdef DEBUG
            assert(i<this->rows);
            assert(i>=0);
            assert(cvPtr2D(this,i,0)==this->data.ptr+i*this->step);
        #endif
        return ((T*)(this->data.ptr+i*this->step))[0];
    }
    T & operator() (const int & i) {
        return ((T*)(this->data.ptr+i*this->step))[0];
    }
    T & operator[] (const int & i) {
        return ((T*)(this->data.ptr+i*this->step))[0];
    }
    const int & size() const { return this->rows; }
    const int & getLength() const { return this->rows; }
    Vector<T> & operator=(const Matrix<T> &m){
        this->resize(m.rows);
		cvCopy(&m, this);
		return (*this);
    }
    Vector<T> & operator=(const Vector<T> &v){
        this->resize(v.size());
		cvCopy(&v, this);
		return (*this);
    }
	Vector<T> & operator=(const T & v){
		cvSet(this, cvScalarAll(v));
		return (*this);
	}
	void resize(const int & i){
		Matrix<T>::resize(i,1);
	}
};
#endif //MATRIX_HH
