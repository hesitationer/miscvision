#ifndef __ARRAY_OPS_HH
#define __ARRAY_OPS_HH

/** Binary operation class */
template <typename array_t, typename func_t> struct Array_BIN_OP {
	const array_t * m_A;
	const array_t * m_B;
	func_t func;
	Array_BIN_OP(const array_t & A, const array_t & B) : 
		m_A(&A), 
		m_B(&B) 
	{
	}
};

// best to understand this by looking at the example below
#define ARRAY_BIN_IMPL(functionname, temp, classname) \
template <typename temp> \
Array_BIN_OP<temp, classname > functionname (const temp & A, const temp & B){ \
    return Array_BIN_OP<temp, classname >(A,B); \
}

#define ARRAY_BIN_IMPL_START(temp, functionname, classname, src1, src2, res) \
template <typename temp> \
struct classname { \
    void operator ()(src1, src2, res) 
	
#define ARRAY_BIN_IMPL_END(temp, functionname, classname) \
}; \
ARRAY_BIN_IMPL(functionname, temp, classname<temp>)

// for defining a C++ overloaded operator C = A op B
#define ARRAY_BIN_OP_START(temp, op, name, src1, src2, res)\
 ARRAY_BIN_IMPL_START(temp, operator op, name##Op, src1, src2, res)
#define ARRAY_BIN_OP_END(temp, op, name)\
 ARRAY_BIN_IMPL_END(temp, operator op, name##Op)

// for defining general functions C=func(A,B) .. where nothing gets copied
#define ARRAY_BIN_FUNC_START(temp, func, src1, src2, res) \
 ARRAY_BIN_IMPL_START(temp, func, func##Func, src1, src2, res)
#define ARRAY_BIN_FUNC_END(temp, func) \
 ARRAY_BIN_IMPL_END(temp, func, func##Func)
	
/*//////////////////// EXAMPLE ///////////////////////////////
//   writing the following code
//////////////////////////////////////////////////////////////
 
ARRAY_BIN_FUNC_START(array_t, add, const array_t & A, const array_t & B, array_t & res)
{
    cvAdd(&A, &B, &res);
}
ARRAY_BIN_FUNC_END(array_t, add) 
 
/////////////////////////////////////////////////////////////
//    generates the following code after expanding macros   /
/////////////////////////////////////////////////////////////
 
 template <typename array_t>
 struct addStruct {
 	void operator ()(const array_t & A, const array_t & B, array_t & res)
	{
		cvAdd(&A, &B, &res);
	}
 };
 template <typename array_t>
 Array_BIN_OP<array_t, addStruct> add (const array_t & A, const array_t & B) {
 	return Array_BIN_OP<array_t, addStruct>(&A, &B);
 }

 this allows one to use the function 'add' as follows
 array_t A, B, C;
 
/////////////////////////////////////////////////////////////
*/

#endif // __ARRAY_OPS_HH
