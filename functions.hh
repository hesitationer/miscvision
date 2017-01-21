#ifndef __FUNCTIONS_HH
#define __FUNCTIONS_HH

#include <cv/array_ops.hh>

#define CV_OP_FUNC_IMPL(op, Func) \
ARRAY_BIN_FUNC_START(pixel_t, Func, const pixel_t & A, const pixel_t & B, pixel_t & res) \
{\
	cv##Func(&A, &B, &res);\
}\
ARRAY_BIN_FUNC_END(pixel_t, Func) \
ARRAY_BIN_OP_START(pixel_t, op, Func, const pixel_t & A, const pixel_t & B, pixel_t & res)\
{\
	cv##Func(&A, &B, &res);\
}\
ARRAY_BIN_OP_END(pixel_t, op, Func) 

// each of these adds a function 
// c = a*b
// c = Mul(a,b)
// etc for +,-...
CV_OP_FUNC_IMPL(+, Add)
CV_OP_FUNC_IMPL(*, Mul)
CV_OP_FUNC_IMPL(-, Sub)
CV_OP_FUNC_IMPL(/, Div)
CV_OP_FUNC_IMPL(&, And)
//CV_OP_FUNC_IMPL(|, Or)
CV_OP_FUNC_IMPL(^, Xor)


#endif // __FUNCTIONS_HH
