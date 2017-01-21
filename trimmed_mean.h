#ifndef __CV_TRIMMED_MEAN_H
#define __CV_TRIMMED_MEAN_H

#include<vector>
#include<algorithm>
#include<cv/matrix.hh> 
template <typename matT, int matS, typename masT, int masS>
CvScalar cvTrimmedMean(const Matrix<matT,matS> & m, int pct, Matrix<masT, masS>* mask){
	char str[256];
	sprintf(str, "Not implemented for ( Matrix<%s,%s>, int, Matrix<%s,%s> * )", STRINGIFY(matT), STRINGIFY(matS) , STRINGIFY(masT) , STRINGIFY(masT) );
	cvError(CV_StsUnsupportedFormat, "cvTrimmedMean", str, __FILE__, __LINE__);
	return cvScalarAll(0);
}

template <typename T>// int matS=1, typename masT=uchar, int masS=1>
CvScalar cvTrimmedMean(const Matrix<T,1> &m, int pct, Matrix<uchar, 1>* mask= NULL){
	std::vector<T> v;
	T * ptr = (T*) m.data.ptr;
	uchar * mptr = NULL;
	int length = m.step*m.rows;
	CvScalar val = cvScalarAll(0);
	
	if(mask) mptr = mask->data.ptr;

	// for each non-zero value, push into v
	if(mask){
		for(int i=0; i<length; i++){
			if(mptr[i] == 0) continue;
			v.push_back( ptr[i] );
		}
	}
	else{
		for(int i=0; i<length; i++){
			v.push_back( ptr[i] );
		}
	}

	// sort values in v
	std::sort(v.begin(), v.end());
	
	// chop last pct percentage
	int nchop = v.size()*pct/100;
	assert(pct<50 && pct>=0);
	v.erase(v.end()-nchop, v.end());
	v.erase(v.begin(), v.begin()+nchop);
	int n=0;
	for(std::vector<T>::iterator it = v.begin();
	    it != v.end();
		++it){
		val.val[0]+=*it;
		n++;
	}
	val.val[0] /= n;
	return val;
}

CvScalar cvTrimmedMean(CvArr * arr, int pct, CvArr * mask=NULL){
	CvScalar result = cvScalarAll(0);
	MATRIX_FROM_CVMAT(arr, mat, MATRIX_FROM_CVMAT(mask,matmask, result=cvTrimmedMean(*mat,pct,matmask)));
	return result;
}

#endif
