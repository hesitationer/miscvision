#ifndef __VECTOR3_HH
#define __VECTOR3_HH

#include <cxcore.h>
#include "cv/cvmat.h"
#include <ostream>
using namespace std;

template <class T> 
struct Vector3{
    T _data[3];

public:
    // Default construction
    Vector3 (const T &x=0, const T &y=0, const T &z=0)
	{
		_data[0] = x;
		_data[1] = y;
		_data[2] = z;
    }
	Vector3 (const CvMat & m){
		(*this) = m;
	}
	Vector3 & operator = (const CvMat & m){
		assert(m.rows>=3);
		this->set(cvmGet(&m,0,0), cvmGet(&m,1,0), cvmGet(&m,2,0));
		return (*this);
	}
    Vector3 (const CvMat * m){
		(*this) = m;
	}
    /*Vector3<T> & operator =(const CvMat * m){
		assert(m!=NULL);
		assert(m->rows>=3);
        this->set(cvmGet(m,0,0), cvmGet(m,1,0), cvmGet(m,2,0));
		return (*this);
    }*/
    Vector3 (const Vector3<T> & v){
        this->set(v(0), v(1), v(2));
    }
	Vector3( const _CvMAT_T_& mat_t ) { *this = CvMAT(mat_t); }
	Vector3( const _CvMAT_INV_& inv_mat ) { *this = CvMAT(inv_mat); }
	Vector3( const _CvMAT_ADD_& mat_add ) { *this = CvMAT(mat_add); }
	Vector3( const _CvMAT_ADD_EX_& mat_add ){ *this = CvMAT(mat_add); }
	Vector3( const _CvMAT_SCALE_& scale_mat ){ *this = CvMAT(scale_mat); }
	Vector3( const _CvMAT_SCALE_SHIFT_& scale_shift_mat ) { *this = CvMAT(scale_shift_mat); }
	Vector3( const _CvMAT_MUL_& mmul ) { *this = CvMAT(mmul); }
	Vector3( const _CvMAT_MUL_ADD_& mmuladd ) { *this = CvMAT(mmuladd); }
	Vector3( const _CvMAT_LOGIC_& mat_logic ) { *this = CvMAT(mat_logic); }
	Vector3( const _CvMAT_UN_LOGIC_& mat_logic ) { *this = CvMAT(mat_logic); }
	Vector3( const _CvMAT_NOT_& not_mat ) { *this = CvMAT(not_mat); }
	Vector3( const _CvMAT_COPY_& mat_copy ) { *this = CvMAT(mat_copy); }
	Vector3( const _CvMAT_CVT_& mat_copy ) { *this = CvMAT(mat_copy); }
	Vector3( const _CvMAT_DOT_OP_& dot_mat ) { *this = CvMAT(dot_mat); }
	Vector3( const _CvMAT_SOLVE_& solve_mat ) { *this = CvMAT(solve_mat); }
	Vector3( const _CvMAT_CMP_& cmp_mat ) { *this = CvMAT(cmp_mat); }
	
    void set(const T &x, const T &y, const T &z){
        _data[0] = x;
        _data[1] = y;
        _data[2] = z;
    }
    T& operator[] (const int & i){
        return _data[i];
    }
    const T& operator[] (const int & i) const {
        return _data[i];
    }
    T& operator() (const int & i){
        return _data[i];
    }
    const T& operator() (const int & i) const {
        return _data[i];
    }
    CvMat * cvMat() const { 
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
        CvMat * m = cvCreateMat(3, 1, type);
        cvmSet(m,0,0, _data[0]);
        cvmSet(m,1,0, _data[1]);
        cvmSet(m,2,0, _data[2]);
        return m;
    }
    friend ostream & operator << (ostream & out, const Vector3 & v){
        return out<<"["<<v(0)<<", "<<v(1)<<", "<<v(2)<<"]";
    }
    Vector3<T> operator * (double x) const {
        return Vector3((T)(_data[0]*x), (T)(_data[1]*x), (T)(_data[2]*x));
    }
    const Vector3<T> & operator *= (double x) {
		_data[0]*=x;
		_data[1]*=x;
		_data[2]*=x;
        return (*this);
    }
    const Vector3<T> & operator *= (Vector3<T> &s) {
		_data[0]*=s._data[0];
		_data[1]*=s._data[1];
		_data[2]*=s._data[2];
        return (*this);
    }
	Vector3<T> operator * (const Vector3<T> & v) const {
		return Vector3(_data[0]*v(0), _data[1]*v(1), _data[2]*v(2));
	}
    Vector3<T> & operator += (const Vector3<T> & v){
		_data[0]+=v[0];
		_data[1]+=v[1];
		_data[2]+=v[2];
        return (*this);
    }
    Vector3<T> & operator += (const T & t){
		_data[0]+=t;
		_data[1]+=t;
		_data[2]+=t;
        return (*this);
    }
	Vector3<T> operator + (const Vector3<T> & v) const{
        return Vector3(_data[0]+v(0), _data[1]+v(1), _data[2]+v(2));
    }
    Vector3<T> operator - (const Vector3<T> & v) const{
        return Vector3(_data[0]-v(0), _data[1]-v(1), _data[2]-v(2));
    }
    Vector3<T> operator - () const{
        return Vector3(-_data[0], -_data[1], -_data[2]);
    }
    const Vector3<T> & operator -= (double x) {
		_data[0]-=x;
		_data[1]-=x;
		_data[2]-=x;
        return (*this);
    }
    const Vector3<T> & operator -= (Vector3<T> &s) {
		_data[0]-=s._data[0];
		_data[1]-=s._data[1];
		_data[2]-=s._data[2];
        return (*this);
    }
    
	Vector3<T> operator / (const Vector3<T> & v){
		return Vector3<T>(_data[0]/v[0], 
						  _data[1]/v[1],
						  _data[2]/v[2]);

	}
	Vector3<T> operator / (const T & x){
		return Vector3<T>(_data[0]/x, 
						  _data[1]/x,
						  _data[2]/x);
	}
	const Vector3<T> & operator /= (const T & x){
		_data[0]/=x;
		_data[1]/=x;
		_data[2]/=x;
		return (*this);
	}
    bool operator == (const T & x){
		return _data[0]==x && _data[1]==x && _data[2]==x;
	}
	Vector3<T> cross(const Vector3<T> & v) const {
        Vector3<T> vNormal;    

        // Calculate the cross product with the non communitive equation
        vNormal[0] = ((_data[1] * v[2]) - (_data[2] * v[1]));
        vNormal[1] = ((_data[2] * v[0]) - (_data[0] * v[2]));
        vNormal[2] = ((_data[0] * v[1]) - (_data[1] * v[0]));

        // Return the cross product
        return vNormal;
    }
    
    const Vector3<T> & operator = (const Vector3<T> & v){
        this->set(v(0), v(1), v(2));   
        return (*this);
    }
    double magnitude() const{
        return sqrt(_data[0]*_data[0]+_data[1]*_data[1]+_data[2]*_data[2]);
    }
    void normalize()
    {
        // Get the magnitude of our normal
        double magnitude = this->magnitude();//Magnitude(vVector);                

        // Now that we have the magnitude, we can divide our vector by that magnitude.
        // That will make our vector a total length of 1.  
        *this = *this/((T)magnitude);
    }
	Vector3<T> normalize() const {
		Vector3<T> x(*this);
		x.normalize();
		return x;
	}
	
    Vector3<T> cross(const Vector3<T> & vVector2)
    {
        Vector3<T> vNormal;

        // Calculate the cross product with the non communitive equation
        vNormal[0] = ((_data[1] * vVector2[2]) - (_data[2] * vVector2[1]));
        vNormal[1]= ((_data[2] * vVector2[0]) - (_data[0] * vVector2[2]));
        vNormal[2] = ((_data[0] * vVector2[1]) - (_data[1] * vVector2[0]));

        // Return the cross product
        return vNormal;
    }
    float dot(const Vector3<T> & v) const{
        return (v[0]*_data[0]+v[1]*_data[1]+v[2]*_data[2]);
    }
	static T dist(const Vector3<T> & a, const Vector3<T> &b){
		return (a-b).magnitude();
	}
};

template <class T> 
struct Vector2{
    T _data[2];

public:
    // Default construction
    Vector2 (const T &x=0, const T &y=0){
        this->set(x,y);
    }
    void set(const T &x, const T &y){
        _data[0] = x;
        _data[1] = y;
    }
    T& operator() (const int & i) {
        return _data[i];
    }
    T& operator[] (const int & i) {
        return _data[i];
    }
    const T& operator() (const int & i) const {
        return _data[i];
    }
    const T& operator[] (const int & i) const {
        return _data[i];
    }
    friend ostream & operator << (ostream & out, const Vector2 & v){
        return out<<"["<<v(0)<<", "<<v(1)<<"]";
    }
};

typedef Vector3<double> dVec3;
typedef Vector3<float>  fVec3;
typedef Vector3<int>    iVec3;
typedef Vector3<uchar>    cVec3;
typedef Vector2<double> dVec2;
typedef Vector2<float>  fVec2;
typedef Vector2<int>    iVec2;
typedef Vector2<uchar>    cVec2;

#endif //__VECTOR3_HH

