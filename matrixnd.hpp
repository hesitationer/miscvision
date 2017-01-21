#ifndef __CV_MATRIXND_HH
#define __CV_MATRIXND_HH

/*M///////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004 Roman Stanchak
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without lmatitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
// permit persons to whom the Software is furnished to do so, subject to the following 
// conditions:
//
// The above copyright notice and this permission notice shall be included in all copies 
// or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
// IN THE SOFTWARE.
//
//M*/

#include <cxcore.h>
#include <iostream>

namespace cv {

// stupid template tricks to get correct depth flag
template <typename T> struct MatrixDepth { static int depth() { assert(0); return 0; }};
template <> struct MatrixDepth<double> {	static int depth() { return CV_64F; }};
template <> struct MatrixDepth<float> {	static int depth() { return CV_32F; }};
template <> struct MatrixDepth<int> {	static int depth() { return CV_32S; }};
template <> struct MatrixDepth<short> {	static int depth() { return CV_16S; }};
template <> struct MatrixDepth<unsigned short> {	static int depth() { return CV_16S; }};
template <> struct MatrixDepth<char> {	static int depth() { return CV_8S; }};
template <> struct MatrixDepth<unsigned char> {	static int depth() { return CV_8U; }};

template <typename T, int ndim, int idim> struct MatrixNDSlice {
	const CvMatND * m_mat;
	uchar *    m_offset;
	
	MatrixNDSlice( const CvMatND * mat, const size_t & idx ):
		m_mat( mat ),
		m_offset( mat->data.ptr + mat->dim[ndim-(idim)].step*idx )
	{
		//	printf("MatrixNDSlice<%d, %d> %p\n", ndim, idim, m_offset);
	}
	MatrixNDSlice( const MatrixNDSlice<T,ndim,idim+1> & slice, const size_t & idx ):
		m_mat( slice.m_mat ),
		m_offset( slice.m_offset + slice.m_mat->dim[ndim-(idim)].step*idx )
	{
			//printf("MatrixNDSlice<%d, %d> %p\n", ndim, idim, m_offset);
	}
	MatrixNDSlice<T,ndim,idim-1> operator[]( const size_t & i ){
		return 	MatrixNDSlice<T,ndim,idim-1>( *this, i );
	}

	MatrixNDSlice<T,ndim,idim-1> operator[]( const size_t & i ) const {
		return 	MatrixNDSlice<T,ndim,idim-1>( *this, i );
	}
};

template <typename T, int ndim> struct MatrixNDSlice<T, ndim, 1>{
	const CvMatND * m_mat;
	uchar *    m_offset;

	// constructors
	MatrixNDSlice( const CvMatND * mat, const size_t & idx ):
		m_mat( mat ),
		m_offset( mat->data.ptr + mat->dim[ndim-1].step*idx )
	{
			//printf("MatrixNDSlice<%d, %d> %p\n", ndim, 2, m_offset);
	}
	MatrixNDSlice( const MatrixNDSlice<T,ndim,2> & slice, const size_t & idx ):
		m_mat( slice.m_mat ),
		m_offset( slice.m_offset + slice.m_mat->dim[ndim-3].step*idx )
	{
			//printf("MatrixNDSlice<%d, %d> %p\n", ndim, 2, m_offset);
	}

	// [] operator functions
	T & operator[]( const size_t & i ){
		return *(T*) (m_offset + m_mat->dim[ndim-1].step*i);
	}
	const T & operator[]( const size_t & i ) const {
		return *(T*) (m_offset + m_mat->dim[ndim-1].step*i);
	}
};

/// C++ wrapper for CvMat structure
template <typename T, int ndim> class MatrixND : public CvMatND {
public:
	
	/** Default Constructor */	
    MatrixND(){
		_setParameters();
    }
	
    /** Constructor.  
	 * @param sizes pointer to an ndim array of dimension sizes 
	 */
    MatrixND(const int * sizes){
		_setParameters();
		this->realloc(sizes);
    }

	MatrixND(int s1, int s2=0, int s3=0, int s4=0){
		int sizes[] = {s1,s2,s3,s4};
		_setParameters();
		this->realloc(sizes);
	}
    
	T * slice( int * index ){
		uchar * ptr=this->data.ptr;
		if(ndim==1) return (T*) ptr;

		for(int i=ndim-2; i>=1; i--){
			ptr+=index[i]*this->dim[i].step;
		}
		ptr+=this->dim[0].step*index[0];
		return ((T*) ptr);
	}
	T * slice(int i=0, int j=0, int k=0){
		return (T*)(this->data.ptr + i*this->dim[0].step + j*this->dim[1].step + k*this->dim[2].step);
	}
	const T & operator() (int i=0, int j=0, int k=0) const{
		return *(T*)(this->data.ptr + i*this->dim[0].step + j*this->dim[1].step + k*this->dim[2].step);
	}
	T & operator() (int i=0, int j=0, int k=0){
		return *(T*)(this->data.ptr + i*this->dim[0].step + j*this->dim[1].step + k*this->dim[2].step);
	}
	MatrixNDSlice<T, ndim, ndim-1>  operator[] (const size_t & i) const{
		return MatrixNDSlice<T, ndim, ndim-1>(this, i);
	}
	MatrixNDSlice<T, ndim, ndim-1>  operator[] (const size_t & i){
		return MatrixNDSlice<T, ndim, ndim-1>(this, i);
	}
	/** set internal data buffer to given pointer */
	void setData(T * data, int * sizes=NULL){
		if(sizes){
			cvInitMatrixNDHeader(this, ndim, sizes, this->depth, (void*) data);
		}
		else{
			cvInitMatrixNDHeader(this, ndim, this->dim, this->depth, (void*) data);
		}
		//this->data.ptr = (char *) data;	
	}
		
	/** allocate or reallocate (if neccessary) the matage buffer.  
	    No memory allocation will take place if w==matage.width and h==matage.height
		@param w width of the reallocd matage
		@param h height of the reallocd matage
	*/
    bool realloc(const int * sizes){
		if(!sizes) return false;
        if(this->data.ptr){
			int i;
			for(i=0; i<ndim; i++){
				if(this->getSize(i)!=sizes[i]) break;
			}

            if(i==ndim){
                return false;
            }
			//printf("release(realloc): %p %d %d\n", this->data.ptr, this->width, this->height);
            cvReleaseData(this);
            data.ptr = NULL;
        }
		_setParameters();
        cvInitMatNDHeader(this, ndim, sizes, MatrixDepth<T>::depth(), NULL );
        cvCreateData(this);

		return true;
    }

	/** allocate or reallocate (if neccessary) the matage buffer.
	  * @param mat CvMat from which the width and height are taken
	  */
    bool realloc (const CvMat & mat){
		int sizes[2] = {mat.width, mat.height};
        return this->realloc(sizes);
    }
	bool realloc (const CvMatND & mat){
		return this->realloc(mat.dim);
	}
    
	/** return the size of the nth dimension */
	inline int getSize(int _dim) const { return this->dim[_dim].size; }
	
	/** return the number of dimensions */
	const int & getND() const { return ndim; }

	~MatrixND(){
		if(this->data.ptr) {
			cvDecRefData( this );
		}
        this->data.ptr = NULL;
	}

private:
	void _setParameters(){
		bzero(this, sizeof(CvMatND));

		this->type = MatrixDepth<T>::depth();
		this->dims = ndim;
		//T mytypevar;
		//cvInitMatNDHeader(this, ndim, sizes, MatrixDepth<T>::depth(), NULL);

	}

};

/* Specialization for case of single dimension Matrix */
template <typename T> class MatrixND<T,1> : public CvMatND {
public:
    const T & operator[] (const size_t & i) const{
		return *(T*)(this->data.ptr + i*this->dim[0].step);
    }
    T & operator[] (const size_t & i){
		return *(T*)(this->data.ptr + i*this->dim[0].step);
    }
};

// Stream operators 
template <class T, int ndim, int idim>
inline std::ostream & operator <<(std::ostream & out, const MatrixNDSlice<T,ndim,idim> &m) {
	int i;
	out<<"[";
	for(i=0; i<m.m_mat->dim[idim].size-1; i++){
		out<<m[i]<<", ";
	}
	out<<m[i];
	out<<"]";
	return out;
}

template <class T, int ndim>
inline std::ostream & operator <<(std::ostream & out, const MatrixND<T,ndim> &m) {
	int i;

	out<<"[";
	for(i=0; i<m.dim[0].size-1; i++){
		out<<m[i]<<",\n";
	}
	out<<m[i];
	out<<"]";
	return out;
}

} // namespace cv
#endif //__CV_MATRIXND_HH

