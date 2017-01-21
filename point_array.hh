#ifndef __POINT_ARRAY_HH
#define __POINT_ARRAY_HH

#include "cv/vector3.hh"
#include "cv/matrix.hh"
#include <cv.h>

template <class T> class PointArray : public Matrix<T> {
public:
	PointArray(const int & len): Matrix<T>(len, 3) {
	}
	
    const Vector3<T> & operator()(const int & i) const{
        return *(Vector3 <T> *)(this->data.ptr + this->step*i); 
    }
    Vector3<T> & operator()(const int & i){
        return *(Vector3 <T> *)(data.ptr + step*i);
    }
	const Vector3<T> & operator[](const int & i) const{
        return *(Vector3 <T> *)(this->data.ptr + this->step*i); 
    }
    Vector3<T> & operator[](const int & i){
        return *(Vector3 <T> *)(data.ptr + step*i);
    }
    T * ptr(const int & i) {
        return (T *)(data.ptr + step*i);
    }
	const int & size() {
		return this->rows;
	}

	/*const PointArray<T> & operator *= (const Matrix<T> & m){
		cvGEMM(this, &m, 1, 0, 1, this);
		return (*this);
	}*/
	const PointArray<T> & operator += (const PointArray<T> &a){
		cvAdd(this,&a,this);
		return (*this);
	}
	const PointArray<T> & operator += (const Vector3<T> & v){
		// cast me to 3 channel matrix
		CvMat ch3;

		int type = this->calcType();
		if(type == CV_32FC1) type = CV_32FC3;
		else if(type == CV_64FC1) type = CV_64FC3;
		else if(type == CV_8UC1) type = CV_8UC3;
		else (type = -1);

		cvInitMatHeader(&ch3, this->size(), 1, type, this->data.ptr);
		cvAddS(&ch3, cvScalar(v[0], v[1], v[2]), &ch3);
		return (*this);
	}
    
	void resize(int size){
        if(size > this->rows){
			cvDecRefData( this );
            int old_size=this->rows*this->step;
            void * ptr = this->data.ptr;
            //this->data.ptr = NULL;
			cvInitMatHeader(this, size, 3, this->calcType());
            cvCreateData(this);
			//memcpy(this->data.ptr, ptr, old_size); 
			// argh! need to free ptr !!!
        }
        else{
            this->rows = size;
        }
    }
	
};

/*
template <class T> class PointArray : public CvMat {
public:
    PointArray(const int & len){
        cvInitMatHeader(this, len, 3, this->calcType());
        cvCreateData(this);
    }
    PointArray(const PointArray & pts){
        cvInitMatHeader(this, pts.rows, 3, pts.type); 
        cvCreateData(this);
        cvCopy(&pts, this);
    }
    PointArray(CvMat * m){
        cvInitMatHeader(this, m->rows, 3, m->type); 
        cvCreateData(this);
        cvCopy(m, this);
    }
    PointArray<T> operator= (const PointArray<T> &pts){
        return PointArray(pts);
    }
    PointArray<T> operator= (CvMat * m){
        return PointArray(m);
    }
    int calcType(){
        int type=CV_32FC1;
        switch (sizeof(T)){
        case 8:
            type = CV_64FC1;
            break;
        case 4:
            type = CV_32FC1;
            break;
        case 1:
            type = CV_8UC1;
            break;
        default:
            assert(0);
            break;
        }
        return type;
    }
    void operator /= (const T & t){
        cvScale(this, this, 1.0/t);
    }
    void operator *= (const T & t){
        cvScale(this, this, t);
    }
    const Vector3<T> & operator()(const int & i) const{
        return *(Vector3 <T> *)(this->data.ptr + this->step*i); 
    }
    Vector3<T> & operator()(const int & i){
        return *(Vector3 <T> *)(data.ptr + step*i);
    }
	const Vector3<T> & operator[](const int & i) const{
        return *(Vector3 <T> *)(this->data.ptr + this->step*i); 
    }
    Vector3<T> & operator[](const int & i){
        return *(Vector3 <T> *)(data.ptr + step*i);
    }
    T * ptr(const int & i) {
        return (T *)(data.ptr + step*i);
    }
	const int & size() {
		return this->rows;
	}
    void resize(int size){
        if(size > this->rows){
            int old_size=this->rows*this->step;
            void * ptr = this->data.ptr;
            this->data.ptr = NULL;
            this->rows = size;
            cvCreateData(this);
            memcpy(this->data.ptr, ptr, old_size); 
			// argh! need to free ptr !!!
        }
        else{
            this->rows = size;
        }
    }
    
    ~PointArray(){
        //free data
        cvDecRefData( this );
    }
};
*/
#endif //__POINT_ARRAY_HH
