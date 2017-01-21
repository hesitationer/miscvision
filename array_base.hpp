#ifndef __ARRAY_BASE_HPP
#define __ARRAY_BASE_HPP

/*M///////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004 Roman Stanchak
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without larritation the rights to use, copy, modify,
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

namespace cv {

// stupid template tricks to get correct depth flag
template <typename T> struct ArrayBaseDepth { static int depth() { assert(0); return 0; }};
template <> struct ArrayBaseDepth<double> {	static int depth() { return CV_64F; }};
template <> struct ArrayBaseDepth<float> {	static int depth() { return CV_32F; }};
template <> struct ArrayBaseDepth<int> {	static int depth() { return CV_32S; }};
template <> struct ArrayBaseDepth<short> {	static int depth() { return CV_16S; }};
template <> struct ArrayBaseDepth<unsigned short> {	static int depth() { return CV_16S; }};
template <> struct ArrayBaseDepth<char> {	static int depth() { return CV_8S; }};
template <> struct ArrayBaseDepth<unsigned char> {	static int depth() { return CV_8U; }};

/// abstract base class for OpenCv array types ... akin to CvArr
template <typename T, int nch=1> class ArrayBase {
public:
	
	/** No Constructor */	

	/** set internal data buffer to given pointer */
	virtual void setData(T * data, int w=-1, int h=-1)=0;
		
	/** return a pointer to the data at row i */
	virtual T* operator[](const int &i) = 0; 

	/** return a constant pointer to the data at row i */
	virtual const T* operator[](const int &i) const = 0;
	
	/** return a constant pointer to the data at row i */
	const T* prow(const int &i) const {
		return (*this)[i];
	}

	/** return a pointer to the data at row i */
	T* rowptr(const int &i) {
		return (*this)[i];
	}

	/** allocate or reallocate (if neccessary) the array buffer.  
	    No memory allocation will take place if w==array.width and h==array.height
		@param w width of the reallocd array
		@param h height of the reallocd array
	*/
    virtual bool realloc(const int & w, const int & h) = 0;

	/** allocate or reallocate (if neccessary) the array buffer.
	  * @param arr CvMat from which the width and height are taken
	  */
    bool realloc (const CvMat & arr){
        return this->realloc(arr.width, arr.height);
    }
	bool realloc (const IplImage & im){
		return this->realloc(im.width, im.height);
	}
    
	/** const individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    virtual inline const T & operator () (const int & x, const int & y, const int &ch=0) const = 0;

	/** individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    virtual inline T & operator () (const int & x, const int & y, const int &ch=0) = 0;

	/** return the width of the array */
	virtual const int & getWidth() const = 0;
	
	/** return the height of the array */
	virtual const int & getHeight() const = 0;
	
	/** alias for cvGetSize */
	inline CvSize getSize() const {
		return cvGetSize(this);
	}

	/** return the depth of the array */
	virtual int getDepth() const = 0;
	
	/** return the number of channels in the array */
	virtual int getNumChannels() const = 0;

		
	/*const ArrayBase<T,nch> & operator= (const CvMat & m){
		assert(arr.data.ptr!=NULL);
		if(arr.data.ptr==NULL) return (*this);
		switch(arr.depth){
			case CV_8U:
				switch(arr.nChannels){
					case 1:
						return this->convert(*(const ArrayBase<unsigned char, 1> *)(&arr));
					case 2:
						return this->convert(*(const ArrayBase<unsigned char, 2> *)(&arr));
					case 3:
						return this->convert(*(const ArrayBase<unsigned char, 3> *)(&arr));
				}
				break;
			case CV_16S:
				switch(arr.nChannels){
					case 1:
						return this->convert(*(const ArrayBase<short, 1> *)(&arr));
					case 2:
						return this->convert(*(const ArrayBase<short, 2> *)(&arr));
					case 3:
						return this->convert(*(const ArrayBase<short, 3> *)(&arr));
				}
				break;
			case CV_32S:
				switch(arr.nChannels){
					case 1:
						return this->convert(*(const ArrayBase<int, 1> *)(&arr));
					case 2:
						return this->convert(*(const ArrayBase<int, 2> *)(&arr));
					case 3:
						return this->convert(*(const ArrayBase<int, 3> *)(&arr));
				}
				break;
			case CV_32F:
				switch(arr.nChannels){
					case 1:
						return this->convert(*(const ArrayBase<float, 1> *)(&arr));
					case 2:
						return this->convert(*(const ArrayBase<float, 2> *)(&arr));
					case 3:
						return this->convert(*(const ArrayBase<float, 3> *)(&arr));
				}
				break;
			case CV_64F:
				switch(arr.nChannels){
					case 1:
						return this->convert(*(const ArrayBase<double, 1> *)(&arr));
					case 2:
						return this->convert(*(const ArrayBase<double, 2> *)(&arr));
					case 3:
						return this->convert(*(const ArrayBase<double, 3> *)(&arr));
				}
				break;
		}
		return (*this);
	}

	template <typename S, int nch2> 
	void eq_operator_same (const ArrayBase<S,nch2> &arr){
		this->realloc(arr);
		cvCopy((CvArr *) &arr, this);
	}
*/	
	/** equals operator */
/*	ArrayBase<T,nch> & operator= (const ArrayBase<T,nch> & arr){
		if(this->data.ptr==NULL && arr.data.ptr==NULL) return (*this);
		assert(arr.data.ptr!=NULL);
		this->realloc(arr.width, arr.height);
		cvCopy(&arr, this);
		return (*this);
	}*/

	/** convert array of possibly different type and color channels */
	/*template <typename S, int nch2> 
	void eq_operator (const ArrayBase<S,nch2> & arr){
		int this_depth=this->getDepth();
		int arr_depth=arr.getDepth();
		int this_nch=nch;
		int arr_nch=nch;
		arr.check();

		if(arr_depth==this_depth){
			if(arr_nch==this_nch){
				// ideally template would take care of this
				this->eq_operator_same( arr );
			}
			else{
				this->eq_operator_nchannels( arr );
			}
		}
		else if(arr_nch==this_nch){
			return this->eq_operator_depth( arr );
		}
		// need to convert twice here, so make sure to do less work
		// by doing type conversion on low color array
		else if(arr_nch > this_nch){
			// convert color first
			ArrayBase<S,nch> dest;
			dest.eq_operator_nchannels(arr);
			// now convert type
			this->eq_operator_depth( dest );
		} 
		else{
			// convert type first
			ArrayBase<T,nch2> dest;
			dest.eq_operator_depth( arr );
			// now convert color 
			this->eq_operator_nchannels( dest );
		} 
	}*/
	
	/** convert array of one type to another -- number of color channels must be the same */
	/*template <typename S, int nch2> 
	void eq_operator_depth (const ArrayBase<S,nch2> & arr){
		this->realloc(arr);
		this->check();
		arr.check();
		CvSize size=arr->getSize();

		assert(nch==nch2); // this should never not happen
		
		for(int j=0; j<size.height; j++){
			const S * s = arr.prow(j);
			T * t = this->prow(j);
			//std::cout<<"t="<<(unsigned long) t<<std::endl;
			//std::cout<<"s="<<(unsigned long) s<<std::endl;
			//ArrayBaseIterator<S> iit = arr.row(j);
			//ArrayBaseIterator<T> tit = this->row(j);
			for(int i=0; i<size.width*nch; i++){
				//tit[i] = (T) iit[i];
				//cvSet2D(this, j,i, cvGet2D(&arr, j, i));
				t[i] = (T) (s[i]);
				//std::cout<<"t[i]="<<(unsigned long) (&(t[i]))<<std::endl;
				//std::cout<<"s[i]="<<(unsigned long) (&(s[i]))<<std::endl;
			}
		}
	}*/

	/** convert color to b/w and vice versa -- just knows about RGB for now */
	/*template <typename S, int nch2> 
	void eq_operator_nchannels (const ArrayBase<S,nch2> & arr){
		int dst_nch = nch;
		int src_nch = nch2;
		
		this->realloc(arr);

		// Many channels -> 1
		if(dst_nch<src_nch && dst_nch==1){
			// I'm lazy, just use the first color channel as the conversion
			// Ideally, this should average the color channels
			cvSplit(&arr, this, NULL, NULL, NULL);
		}
		// One channel -> many
		else if(src_nch < dst_nch && src_nch==1){
			CvArr * channels[4];
			for(int i=0; i<4; i++){
				if(i<dst_nch){
					channels[i]=(CvArr*)&arr;
				}
				else{
					channels[i]=NULL;
				}
			}
			cvMerge(channels[0], channels[1], channels[2], channels[3], this);
		}
		// Otherwise action is undefined
		else{
			assert(0);
		}
	}
	*/
};

} // namespace cv


#endif //__ARRAY_BASE_HPP
