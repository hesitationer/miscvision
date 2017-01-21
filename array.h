#ifndef __CV_ARRAY_HH
#define __CV_ARRAY_HH

#include <iostream>

namespace cv {
template <typename T, int length> struct Array;

// inductive case
#define array_operation_impl_n(func, op)                \
static inline void func(const Array<T,length> & x,     \
		const Array<T,length> & y,                    \
		Array<T, length> &z){                         \
	z._data[index] = x._data[index] op y._data[index]; \
	array_operation<T,length,index-1>::func(x,y,z);    \
}                                                     \
static inline void func(const Array<T,length> & x,     \
		Array<T, length> &z){                         \
	z._data[index] = z._data[index] op x._data[index]; \
	array_operation<T,length,index-1>::func(x,z);    \
}                                                      \
static inline void func(const Array<T,length> & x,     \
		const T & y,                                  \
		Array<T, length> &z){                         \
	z._data[index] = x._data[index] op y;              \
	array_operation<T,length,index-1>::func(x,y,z);    \
}                                                     \
static inline void func(const T & x,                   \
		Array<T, length> &z){                         \
	z._data[index] = z._data[index] op x._data[index];                 \
	array_operation<T,length,index-1>::func(x,z);      \
}

// for logical operations
#define array_logic_operation_impl_n(func, op)               \
static inline bool func(const Array<T,length> & x, const Array<T,length> & y) { \
	return (x[index] op y[index]) && array_operation<T,length,index-1>::func(x,y); \
}\
static inline bool func(const Array<T,length> & x, const T & y) { \
	return (x[index] op y) && array_operation<T,length,index-1>::func(x,y); \
}

// base case
#define array_operation_impl_0(func, op)                \
static inline void func(const Array<T,length> & x,     \
		const Array<T,length> & y,                    \
		Array<T, length> &z){                         \
	z._data[0] = x._data[0] op y._data[0]; \
}                                                     \
static inline void func(const Array<T,length> & x,     \
		Array<T, length> &z){                         \
	z._data[0] = z._data[0] op x._data[0]; \
}                                                      \
static inline void func(const Array<T,length> & x,     \
		const T & y,                                  \
		Array<T, length> &z){                         \
	z._data[0] = x._data[0] op y;              \
}                                                     \
static inline void func(const T & x,                   \
		Array<T, length> &z){                         \
	z._data[0] = z._data[0] op x._data[0];                 \
}

#define array_logic_operation_impl_0(func, op)               \
static inline bool func(const Array<T,length> & x, const Array<T,length> & y) { \
	return x[0] op y[0]; \
}\
static inline bool func(const Array<T,length> & x, const T & y) { \
	return x[0] op y; \
}

// these are defines for operators on arrays, for instance,
// array_operator_impl(add, +, +=) produces
// array_t operator+(array_t)
// array_t operator+(T)
// array_t operator+=(array_t)
// array_t operator+=(array_t)
#define array_operator_impl(func, op, ope) \
	array_t operator op (const array_t & x) const{ \
		array_t z; \
		array_operation<T,length,length-1>::func((*this),x,z); \
		return z;\
	}\
	array_t operator op (const T & x) const{ \
		array_t z; \
		array_operation<T,length,length-1>::func((*this),x,z); \
		return z;\
	}\
	array_t operator ope (const array_t & x) { \
		array_operation<T,length,length-1>::func((*this),x,(*this)); \
		return (*this);\
	}\
	array_t operator ope (const T & x) { \
		array_operation<T,length,length-1>::func((*this),x,(*this)); \
		return (*this);\
	}

#define array_logic_operator_impl(func, op) \
	bool operator op (const array_t & x) const{ \
		return array_operation<T,length,length-1>::func((*this), x);\
	}\
	bool operator op (const T & x) const { \
		return array_operation<T,length,length-1>::func((*this), x); \
	}

// inductive case
template <typename T, int length, int index> 
class  array_operation {
	public:
	array_operation_impl_n(add, +)
	array_operation_impl_n(sub, -)
	array_operation_impl_n(mul, *)
	array_operation_impl_n(div, /)
	array_logic_operation_impl_n(gt, >)
	array_logic_operation_impl_n(lt, <)
	array_logic_operation_impl_n(eqe, ==)
	array_logic_operation_impl_n(ne, !=)
	array_logic_operation_impl_n(gte, >=)
	array_logic_operation_impl_n(lte, <=)
	static inline std::ostream& print(std::ostream & o, const Array<T,length> & x){
		o<<x._data[length-index-1]<<", ";
		return array_operation<T,length,index-1>::print(o,x);
	}
	static inline void eq(const Array<T,length> & x,
	                      Array<T,length> &z){
		z._data[index] = x._data[index];
		array_operation<T,length,index-1>::eq(x,z);
	}
	static inline void eq(const T & x,
	                      Array<T,length> &z){
		z._data[index] = x;
		array_operation<T,length,index-1>::eq(x,z);
	}
	static inline void sum(const Array<T,length> & x, double & res){
		res += x._data[index];
		array_operation<T,length,index-1>::sum(x,res);
	}
};

// base case
template <typename T, int length>
class array_operation<T, length, 0> {
public:
	array_operation_impl_0(add, +)
	array_operation_impl_0(sub, -)
	array_operation_impl_0(mul, *)
	array_operation_impl_0(div, /)
	array_logic_operation_impl_0(gt, >)
	array_logic_operation_impl_0(lt, <)
	array_logic_operation_impl_0(eqe, ==)
	array_logic_operation_impl_0(ne, !=)
	array_logic_operation_impl_0(gte, >=)
	array_logic_operation_impl_0(lte, <=)
	static inline std::ostream& print(std::ostream & o, const Array<T,length> & x){
		return o<<x._data[length-1];
	}
	static inline void eq(const Array<T,length> & x,
	                      Array<T,length> &z){
		z._data[0] = x._data[0];
	}
	static inline void eq(const T & x,
	                      Array<T,length> &z){
		z._data[0] = x;
	}
	static inline void sum(const Array<T,length> & x, double & res){
		res += x._data[0];
	}
};

template <typename T, int length> struct Array {
	typedef Array<T, length> array_t;
	
	T _data[length];

	Array(){
	}
	Array (T & x1){
		_data[0] = x1;
	}
	Array (T & x1, T &x2){
		_data[0] = x1;
		_data[1] = x2;
	}
	Array (T & x1, T &x2, T &x3){
		_data[0] = x1;
		_data[1] = x2;
		_data[2] = x3;
	}
	Array (T & x1, T &x2, T &x3, T &x4){
		_data[0] = x1;
		_data[1] = x2;
		_data[2] = x3;
		_data[3] = x4;
	}

	// constructor, assumes data is of length 'length'
	Array(T* p){
		memcpy(_data, p, sizeof(T)*length);
	}
	Array(const Array & a){
		array_operation<T,length,length-1>::eq(a, (*this));
	}
	array_t & operator= (const array_t & a){
		array_operation<T,length,length-1>::eq(a, (*this));
		return (*this);
	}
	const T* get_ptr() const {
		return _data;
	}
	/*array_t & operator= (const float &x){
		array_operation<T,length,length-1>::eq(x, (*this));
		return (*this):
	}*/
	array_t & operator= (const T & x){
		array_operation<T,length,length-1>::eq(x, (*this));
		return (*this);
	}
	inline T & operator [] (const int & i) {
		return _data[i];
	}
	const inline T & operator [] (const int & i) const{
		return _data[i];
	}
	array_operator_impl(add,+,+=)
	array_operator_impl(sub,-,-=)
	array_operator_impl(mul,*,*=)
	array_operator_impl(div,/,/=)
	
	array_logic_operator_impl(gt, >)
	array_logic_operator_impl(lt, <)
	array_logic_operator_impl(eqe, ==)
	array_logic_operator_impl(ne, !=)
	array_logic_operator_impl(gte, >=)
	array_logic_operator_impl(lte, <=)
	
	friend std::ostream & operator << (std::ostream & out, const array_t & a){
		out<<"[";
		array_operation<T,length,length-1>::print(out,a);
		return out<<"]";
	}
	void normalize() {
		(*this)/=this->sum();
	}
	double sum() const {
		double result=0;
		array_operation<T,length,length-1>::sum(*this,result);
		return result;
	}
	int size() const {
		return length;
	}

};

}

#endif //__CV_ARRAY_HH
