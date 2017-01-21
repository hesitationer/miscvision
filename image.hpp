#ifndef __CV_IMAGE_HH
#define __CV_IMAGE_HH

/*M///////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2004 Roman Stanchak
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without limitation the rights to use, copy, modify,
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
#include <cv/array_base.hpp>

namespace cv {

// stupid template tricks to get correct IPL depth flag
template <typename T> struct ImageDepth { static int depth() { assert(0); return 0; }};
template <> struct ImageDepth<double> {	static int depth() { return IPL_DEPTH_64F; }};
template <> struct ImageDepth<float> {	static int depth() { return IPL_DEPTH_32F; }};
template <> struct ImageDepth<int> {	static int depth() { return IPL_DEPTH_32S; }};
template <> struct ImageDepth<short> {	static int depth() { return IPL_DEPTH_16S; }};
template <> struct ImageDepth<unsigned short> {	static int depth() { return IPL_DEPTH_16S; }};
template <> struct ImageDepth<char> {	static int depth() { return IPL_DEPTH_8S; }};
template <> struct ImageDepth<unsigned char> {	static int depth() { return IPL_DEPTH_8U; }};

/// C++ wrapper for IplImage structure
template <typename T, int nch=1> class Image : public IplImage { //public ArrayBase<T, nch> {
public:
	
	/** Default Constructor */	
    Image(){
		_setParameters();
    }
	
    /** Constructor.  
	 * @param w width of new image
	 * @param h height of new image
	 */
    Image(int w, int h){
		_setParameters();
		this->realloc(w,h);
    }
    
    /** Copy constructor */
	Image(const Image<T,nch> & im){
        _setParameters();
		this->realloc(im);
        this->operator=(im);
    }

	/** Copy Constructor for IplImage */
	Image(const IplImage & im){
		//assert("cv::Image(const IplImage & im) is broken"==0);
		_setParameters();
		this->realloc(im);
		this->operator=(im);
	}
	
	/** set Region of Interest */
	void setImageROI(const CvRect &r) {
		cvSetImageROI(this, r);
	}

	void resetImageROI() {
		cvResetImageROI(this);
	}
	
	/** set internal image buffer to given pointer */
	void setData(T * data, int w=-1, int h=-1){
		if(w!=-1 && h!=-1){
			cvInitImageHeader(this, cvSize(w,h), this->depth, this->nChannels);
		}
		this->imageData = (char *) data;	
	}
		
	/** return a pointer to the data at row i */
	T* operator[](const int &i) {
		return (T*) (this->imageData+(i*this->widthStep));
	}

	/** return a constant pointer to the data at row i */
	const T* operator[](const int &i) const {
		return (T*) (this->imageData+(i*this->widthStep));
	}
	
	/** allocate or reallocate (if neccessary) the image buffer.  
	    No memory allocation will take place if w==image.width and h==image.height
		@param w width of the reallocd image
		@param h height of the reallocd image
	*/
    bool realloc(const int & w, const int & h){
		if(w<=0 || h<=0) return false;
        if(this->imageData){
            if(this->width==w && 
               this->height==h){ 
                return false;
            }
			//printf("release(realloc): %p %d %d\n", this->imageData, this->width, this->height);
            cvReleaseData(this);
            imageData = NULL;
        }
		_setParameters();
        cvInitImageHeader(this, cvSize(w,h), this->depth, this->nChannels);
		printf("%d %d %d\n", sizeof(*this), sizeof(IplImage), this->nSize);
		assert(this!=0);
		assert(this->nSize==sizeof(IplImage));
		assert(CV_IS_IMAGE_HDR( this ) );
        cvCreateData(this);
		//printf("create(realloc): %p %d %d\n", this->imageData, this->width, this->height);
		return true;
    }

	/** const individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    inline const T & operator () (const int & x, const int & y, const int &ch=0) const {
        return *(T *)(this->imageData + y*this->widthStep + (x*this->nChannels + ch)*sizeof(T));    
        //return ((T *) this->imageData)[(y*this->width + x)*this->nChannels + ch];    
    }

	/** individual pixel access.  It is preferable to grab a pointer to the whole row 
	    and iterate over that */
    inline T & operator () (const int & x, const int & y, const int &ch=0){
        return *(T *)(this->imageData + y*this->widthStep + (x*this->nChannels + ch)*sizeof(T));    
        //return ((T *) this->imageData)[(y*this->width + x)*this->nChannels + ch];    
    }

	/** return the width of the image */
	const int & getWidth() const { return this->width; }
	
	/** return the height of the image */
	const int & getHeight() const { return this->height; }
	
	/** return the depth of the image */
	int getDepth() const { return this->depth; }
	
	/** return the number of channels in the image */
	int getNumChannels() const { return this->nChannels; }

	
	virtual ~Image(){
		if(this->imageData) {
			cvReleaseData(this);
		}
        this->imageData = NULL;
	}

private:
	void _setParameters(){
		//T mytypevar;
		cvInitImageHeader(this, cvSize(0,0), ImageDepth<T>::depth(), nch, IPL_ORIGIN_TL, 4);
		this->imageData=NULL;
	}
};

} // namespace cv

#endif //__CV_IMAGE_HH
