#ifndef __IMAGE_CVARR_ADAPTER_HPP
#define __IMAGE_CVARR_ADAPTER_HPP

#include <cxerror.h>
#include <cv/image.h>

// Often want to convert function of type 
// void func( cv::Image & im )
// to the more generic form
// void func( CvArr * arr)
//
// The classes below will allow you to so as follows:
//
// given void func( cv::Image & im)
// define the struct
//
// struct func_cvarr_adapter {
//    int operator() ( cv::Image & im ){
// 			func( im );
//    }
// };
//
// then func( CvArr * arr) should look as follows
//
// void func( CvArr * arr ) {
// 	  IplImageAdapter1< func_cvarr_adapter > adapter_func;
// 	  adapter_func( arr );
// }
//
// for passing other arguments, the best current method is adding member
// variables to func_cvarr_adapter struct
//

#define IPL_FUNC_SWITCH_ENTRY( full_func_call, imname, ipl_arg, data_t, nch ) \
	{ \
		cv::Image<data_t, nch> & imname = *(cv::Image<data_t, nch> *) ipl_arg; \
		return full_func_call; \
	} 

#define IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, data_t ) \
	switch( ipl_arg->nChannels ){ \
		case 1: \
			IPL_FUNC_SWITCH_ENTRY( full_func_call, imname, ipl_arg, data_t, 1)\
		case 2: \
			IPL_FUNC_SWITCH_ENTRY( full_func_call, imname, ipl_arg, data_t, 2)\
		case 3: \
			IPL_FUNC_SWITCH_ENTRY( full_func_call, imname, ipl_arg, data_t, 3)\
		case 4: \
			IPL_FUNC_SWITCH_ENTRY( full_func_call, imname, ipl_arg, data_t, 4)\
		default: \
			CV_ERROR(CV_BadNumChannels, ""); \
	}

#define IPL_FUNC_SWITCH_VERBOSE( full_func_call, cvarr, imname) \
	IplImage ipl_header; \
	IplImage * ipl_arg = cvGetImage( cvarr, &ipl_header); \
	switch( ipl_arg->depth ){ \
		case IPL_DEPTH_8S: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, char ) \
			break; \
		case IPL_DEPTH_8U: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, uchar ) \
			break; \
		case IPL_DEPTH_16S: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, short ) \
			break; \
		case IPL_DEPTH_16U: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, ushort ) \
			break; \
		case IPL_DEPTH_32S: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, int ) \
			break; \
		case IPL_DEPTH_32F: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, float ) \
			break; \
		case IPL_DEPTH_64F: \
			IPL_FUNC_SWITCH_NCH( full_func_call, imname, ipl_arg, double ) \
			break; \
		default: \
			CV_ERROR(CV_BadDepth, "");\
	}

#define IPL_FUNC_SWITCH( func, cvarr ) \
	IPL_FUNC_SWITCH_VERBOSE( func( ipl_##cvarr ), cvarr, ipl_##cvarr )

template <typename data_t, typename func_t>
struct CvMatFuncSwitch {
	static int apply( CvMat * mat, func_t & func ){
		CV_FUNCNAME("CvMatFuncSwitch::apply");
		__BEGIN__
		switch(CV_MAT_CN(mat->type)) {
			case 1:
				return func( * (Matrix<data_t, 1> *) mat );
			case 2:
				return func( * (Matrix<data_t, 2> *) mat );
			case 3:
				return func( * (Matrix<data_t, 3> *) mat );
			case 4:
				return func( * (Matrix<data_t, 4> *) mat );
			default: 
				CV_ERROR(CV_BadNumChannels, "");
		}
		__END__
		return 0;
	}
	static int apply( CvMat * mat1, CvMat * mat2, func_t & func ){
		CV_FUNCNAME("CvMatFuncSwitch::apply");
		__BEGIN__
		switch(CV_MAT_CN(mat1->type)) {
			case 1:
				return func( * (Matrix<data_t, 1> *) mat1, * Matrix<data_t, 1>::safe_cast(mat2) );
			case 2:
				return func( * (Matrix<data_t, 2> *) mat1, * Matrix<data_t, 2>::safe_cast(mat2) );
			case 3:
				return func( * (Matrix<data_t, 3> *) mat1, * Matrix<data_t, 3>::safe_cast(mat2) );
			case 4:
				return func( * (Matrix<data_t, 4> *) mat1, * Matrix<data_t, 4>::safe_cast(mat2) );
			default: 
				CV_ERROR(CV_BadNumChannels, "");
		}
		__END__
		return 0;
	}
};

template <typename data_t, typename func_t>
struct IplFuncSwitch {
	static int apply( IplImage * im, func_t & func ){
		CV_FUNCNAME("IplFuncSwitch::apply");
		__BEGIN__
		switch(im->nChannels) {
			case 1:
				return func( * (cv::Image<data_t, 1> *) im );
			case 2:
				return func( * (cv::Image<data_t, 2> *) im );
			case 3:
				return func( * (cv::Image<data_t, 3> *) im );
			case 4:
				return func( * (cv::Image<data_t, 4> *) im );
			default: 
				CV_ERROR(CV_BadNumChannels, "");
		}
		__END__
		return 0;
	}
};

template <typename func_t>
struct CvMatAdapter1 {
	func_t func;
	CvMatAdapter1()
	{}

	CvMatAdapter1( const func_t & f) :
		func( f )
	{}

	int operator () (CvArr * cvarr1, CvArr * cvarr2){
		int ret;
		CV_FUNCNAME("CvMatAdapter1");
		__BEGIN__
		CvMat mat1_header; 
		CvMat * mat1_arg = cvGetMat( cvarr1, &mat1_header);
		CvMat mat2_header; 
		CvMat * mat2_arg = cvGetMat( cvarr2, &mat2_header);

		switch( CV_MAT_DEPTH( mat1_arg->type ) ){
		case CV_8S: 
			ret = CvMatFuncSwitch<char, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_8U: 
			ret = CvMatFuncSwitch<uchar, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_16S: 
			ret = CvMatFuncSwitch<short, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_16U: 
			ret = CvMatFuncSwitch<ushort, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_32S: 
			ret = CvMatFuncSwitch<int, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_32F: 
			ret = CvMatFuncSwitch<float, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		case CV_64F: 
			ret = CvMatFuncSwitch<double, func_t>::apply( mat1_arg, mat2_arg, func );
			break;
		default: 
			 CV_ERROR(CV_BadDepth, "");
    	}

		__END__
		return ret;
	}
	int operator () (CvArr * cvarr){
		int ret;
		CV_FUNCNAME("CvMatAdapter1");
		__BEGIN__
		CvMat mat_header; 
		CvMat * mat_arg = cvGetMat( cvarr, &mat_header); 
		switch( CV_MAT_DEPTH( mat_arg->type ) ){
		case CV_8S: 
			ret = CvMatFuncSwitch<char, func_t>::apply( mat_arg, func );
			break;
		case CV_8U: 
			ret = CvMatFuncSwitch<uchar, func_t>::apply( mat_arg, func );
			break;
		case CV_16S: 
			ret = CvMatFuncSwitch<short, func_t>::apply( mat_arg, func );
			break;
		case CV_16U: 
			ret = CvMatFuncSwitch<ushort, func_t>::apply( mat_arg, func );
			break;
		case CV_32S: 
			ret = CvMatFuncSwitch<int, func_t>::apply( mat_arg, func );
			break;
		case CV_32F: 
			ret = CvMatFuncSwitch<float, func_t>::apply( mat_arg, func );
			break;
		case CV_64F: 
			ret = CvMatFuncSwitch<double, func_t>::apply( mat_arg, func );
			break;
		default: 
			 CV_ERROR(CV_BadDepth, "");
    	}

		__END__
		return ret;
	}
};

template <typename func_t>
struct IplImageAdapter1 {
	func_t func;
	IplImageAdapter1()
	{}

	IplImageAdapter1( const func_t & f) :
		func( f )
	{}

	int operator () (CvArr * cvarr){
		int ret;
		CV_FUNCNAME("IplImageAdapter1");
		__BEGIN__
		IplImage ipl_header; 
		IplImage * ipl_arg = cvGetImage( cvarr, &ipl_header); 
		switch( ipl_arg->depth ){ 
		case IPL_DEPTH_8S: 
			ret = IplFuncSwitch<char, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_8U: 
			ret = IplFuncSwitch<uchar, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_16S: 
			ret = IplFuncSwitch<short, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_16U: 
			ret = IplFuncSwitch<ushort, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_32S: 
			ret = IplFuncSwitch<int, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_32F: 
			ret = IplFuncSwitch<float, func_t>::apply( ipl_arg, func );
			break;
		case IPL_DEPTH_64F: 
			ret = IplFuncSwitch<double, func_t>::apply( ipl_arg, func );
			break;
		default: 
			 CV_ERROR(CV_BadDepth, "");
    	}

		__END__
		return ret;
	}
};

template <typename func_t>
struct IplImageAdapter2 {
	func_t func;
	IplImageAdapter2()
	{}

	IplImageAdapter2( const func_t & f) :
		func( f )
	{}

	int operator () (CvArr * im1, CvArr * im2){
		CV_FUNCNAME("IplImageAdapter2");
		__BEGIN__
		IPL_FUNC_SWITCH_VERBOSE( (*this)(cvim1, im2), im1, cvim1 );
		__END__
		return 0;
	}

	template <typename T, int nch>
	int operator () (cv::Image<T, nch> & im1, CvArr * im2){
		CV_FUNCNAME("IplImageAdapter2");
		__BEGIN__
		IPL_FUNC_SWITCH_VERBOSE( func(im1, cvim2), im2, cvim2 );
		__END__
		return 0;
	}
};

// The simple case -- function must be of form  
//  	template <typename T, int nch>
// 		int func(cv::Image<T,nch> & im)
#define CVARR_TO_IMAGE_FUNC1( image_func )\
struct image_func##_cvarr_adapter { \
	template<typename T, int nch> \
	int operator() ( cv::Image<T, nch> & im ) { \
		return image_func( im ); \
	} \
}; \
int image_func( CvArr * arr ){ \
	IplImageAdapter1< image_func##_cvarr_adapter > func_adapter; \
	return func_adapter( arr ); \
}

// same thing for Matrix<T,nch>
#define CVARR_TO_MAT_FUNC_UNARY( mat_func )\
struct mat_func##_cvarr_adapter { \
	template<typename T, int nch> \
	int operator() ( Matrix<T, nch> & mat ) { \
		return mat_func( mat ); \
	} \
}; \
int mat_func( CvArr * arr ){ \
	CvMatAdapter1< mat_func##_cvarr_adapter > func_adapter; \
	return func_adapter( arr ); \
}

#define CVARR_TO_MAT_FUNC_BINARY( mat_func )\
struct mat_func##_cvarr_adapter { \
	template<typename T, int nch> \
	int operator() ( Matrix<T, nch> & mat, Matrix<T, nch> & mat2 ) { \
		return mat_func( mat, mat2 ); \
	} \
}; \
int mat_func( CvArr * arr , CvArr * arr2 ){ \
	CvMatAdapter1< mat_func##_cvarr_adapter > func_adapter; \
	return func_adapter( arr, arr2 ); \
}

// Two arguments -- function of form
//  	template <typename T1, int nch1, typename T2, int nch2>
// 		int func(cv::Image<T,nch> & im1, cv::Image<T, nch> & im2)
#define CVARR_TO_IMAGE_FUNC2( image_func )\
struct image_func##_cvarr_adapter { \
	template<typename T1, int nch1, typename T2, int nch2> \
	int operator() ( cv::Image<T1, nch1> & im1, cv::Image<T2, nch2> & im2) { \
		return image_func( im1, im2 ); \
	} \
}; \
int image_func( CvArr * arr1, CvArr * arr2 ){ \
	IplImageAdapter2< image_func##_cvarr_adapter > func_adapter; \
	return func_adapter( arr1, arr2 ); \
}

#endif // __IMAGE_CVARR_ADAPTER_HPP
