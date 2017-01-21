MAIN: {
#include <cxcore.h>	
print<<END_OF_TEXT;
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
END_OF_TEXT
	gen_array_binary_operator_impl("+", "add", "cvAdd");
	gen_array_binary_operator_impl("-", "sub", "cvSub");
	gen_array_binary_operator_impl("*", "mul", "cvMul");
	gen_array_binary_operator_impl("/", "div", "cvDiv");
	gen_array_binary_operator_impl("^", "xor", "cvXor");
	gen_array_binary_operator_impl("|", "or", "cvOr");
	gen_array_binary_operator_impl("7", "and", "cvAnd");
}
sub gen_array_binary_function {
	my ($function, $type, $class) = @_;
	print "template <typename $type>\n";
	print "Array_BIN_OP< $type, $class > $function(const $type & A, const $type &B){\n";
	print "\treturn Array_BIN_OP< $type, $class >(A,B);\n";
	print "}\n";
}

sub gen_array_binary_function_struct_header {
	my ($temp, $functionname, $classname, $src1, $src2, $res) = @_;
	print "template <typename $temp>\n";
	print "struct $classname { \n";
    print "\tvoid operator ()($src1, $src2, $res)\n";
}

sub gen_array_binary_function_struct_footer {
	my ($temp, $functionname, $classname) = @_;
	print "};\n";
	gen_array_binary_function($functionname, $temp, "$classname< $temp >");
}
sub gen_array_binary_operator_impl {
	my ($operator, $functionname, $functiontocall) = @_;
	gen_array_binary_function_struct_header("array_t", $functionname, $functiontocall."Func", "const array_t & A", "const array_t & B", "array_t & res");
	print "\t{\n\t\t$functiontocall(&A, &B, &res);\n\t}\n";	
	print "};\n";
	
	gen_array_binary_function("operator$operator", "array_t", $functiontocall."Func<array_t>");
	gen_array_binary_function($functionname, "array_t", $functiontocall."Func<array_t>");
	
}
sub gen_array_binary_function_impl {
	my ($functionname, $functiontocall) = @_;
	gen_array_binary_function_struct_header("array_t", $functionname, $functiontocall."Func", "const array_t & A", "const array_t & B", "array_t & res");
	print "\t{\n\t\t$functiontocall(&A, &B, &res);\n\t}\n";	
	gen_array_binary_function_struct_footer("array_t", $functionname, $functiontocall."Func");
}
