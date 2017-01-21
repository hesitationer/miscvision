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

// enable compilation of C++ Image functions
#ifndef USE_CPP_IMAGE
#define USE_CPP_IMAGE
#endif

#include <fstream>
#include <iostream>
#include <cxcore.h>
#include <cv.h>
#include <highgui.h>
//#include <cv/matrix.hh>

template <typename T, int nch> class Matrix;

namespace cv {

template <typename pixel_t, int nch> struct Image_OP;
template <typename pixel_t, int nch> class Image;

// stupid template tricks to get correct IPL depth flag
template <typename T> struct ImageDepth { static int depth() { assert(0); return 0; }};
template <> struct ImageDepth<double> {	static int depth() { return IPL_DEPTH_64F; }};
template <> struct ImageDepth<float> {	static int depth() { return IPL_DEPTH_32F; }};
template <> struct ImageDepth<int> {	static int depth() { return IPL_DEPTH_32S; }};
template <> struct ImageDepth<short> {	static int depth() { return IPL_DEPTH_16S; }};
template <> struct ImageDepth<unsigned short> {	static int depth() { return IPL_DEPTH_16S; }};
template <> struct ImageDepth<char> {	static int depth() { return IPL_DEPTH_8S; }};
template <> struct ImageDepth<unsigned char> {	static int depth() { return IPL_DEPTH_8U; }};

template <typename T, int nch=1> class ImageHeader : public IplImage {
public:
	ImageHeader(){
		_setParameters();
	}
	ImageHeader(int w, int h){
		_setParameters();
	}
	ImageHeader(CvArr * arr){
		if(CV_IS_IMAGE(arr)){
			memcpy(this, arr, sizeof(IplImage));
		}
		else{
			cvGetImage( arr, this );
		}
	}
	operator const Image<T, nch> & () const{
		return *(Image<T,nch> *) this;
	}
	operator Image<T, nch> & (){
		return *(Image<T,nch> *) this;
	}


protected:
	void _setParameters(){
		//T mytypevar;
		cvInitImageHeader(this, cvSize(0,0), ImageDepth<T>::depth(), nch, IPL_ORIGIN_TL, 4);
		this->imageData=NULL;
	}
};

/// C++ wrapper for IplImage structure
template <typename T, int nch=1> class Image : public ImageHeader<T, nch>{
public:
	
	/** Default Constructor */	
    Image(){
		ImageHeader<T,nch>::_setParameters();
    }
	
    /** Constructor.  
	 * @param w width of new image
	 * @param h height of new image
	 */
    Image(int w, int h){
		ImageHeader<T,nch>::_setParameters();
		this->realloc(w,h);
    }
    
    /** Copy constructor */
	Image(const Image<T,nch> & im){
		//assert("copy constructor is broken"==0);
        ImageHeader<T,nch>::_setParameters();
		this->realloc(im);
        this->operator=(im);
    }
	/** Copy constructor for different template type */
	/*template <typename S>
	Image(const Image<S> & im){
		_setParameters();
		this->realloc(im);
		(*this) = im;
	}*/

	/** Copy Constructor for IplImage */
	Image(const IplImage & im){
		//assert("cv::Image(const IplImage & im) is broken"==0);
		ImageHeader<T,nch>::_setParameters();
		this->realloc(im);
		this->operator=(im);
	}
	
	/** explicit image type conversion to char with scaling */
	Image<unsigned char> imagesc(const char * name) const{
		double min,max;
		cvMinMaxLoc(this, &min, &max);
		return this->imagesc(name, min, max);
		
	}	
	Image<unsigned char> imagesc(const char * name, float min, float max) const{
        Image<unsigned char> im;
		im.realloc(this);
		if(min==max){
			cvZero(&im);
		}
		else{
			cvScale(this, &im, 255.0/(max-min), 255*(-min)/(max-min));
		}
		im.show(name);
		return im;
	}
	/** set Region of Interest */
	void setImageROI(const CvRect &r) {
		cvSetImageROI(this, r);
	}
	void resetImageROI() {
		cvResetImageROI(this);
	}
	CvRect getImageROI() const {
		return cvGetImageROI(this);
	}
	
	/** Extract part of image */
	Image<T,nch> getSubImage(CvPoint * pts) const{
		int maxX = MAX(pts[0].x, pts[1].x);
		int minX = MIN(pts[0].x, pts[1].x);
		int maxY = MAX(pts[0].y, pts[1].y);
		int minY = MIN(pts[0].y, pts[1].y);
		
		return this->getSubImage(cvRect(minX,minY,maxX-minX, maxY-minY));
	}
	void getSubImage(CvRect r, Image<T,nch> & dest) const{
		CvMat m,me;
		 cvGetMat( this, &m );
		 cvGetMat( this, &me );
		 r.x = MAX(r.x, 0);
		 r.y = MAX(r.y, 0);
		 r.width = MIN(this->width-r.x, r.width);
		 r.height = MIN(this->height-r.y, r.height);
		 assert(r.width>=0);
		 assert(r.height>=0);
		 cvGetSubRect(&me, &m, r);
		 dest.realloc(r.width,r.height);
		 cvCopy(&m, &dest); 
	}
	Image<T,nch> getSubImage(const CvRect & r) const{
		Image<T,nch> im(r.width, r.height);
		this->getSubImage(r, im);
		return im;
	}

	/** set internal image buffer to given pointer */
	void setData(T * data, int w=-1, int h=-1){
		if(w!=-1 && h!=-1){
			cvInitImageHeader(this, cvSize(w,h), this->depth, this->nChannels);
		}
		this->imageData = (char *) data;	
		this->imageDataOrigin = (char *) data;	
	}
		
	/** scale the image intensity linearly such that inmin maps to outmin and 
	 *  inmax maps to outmax */
	void convert_scale(double inmin, double inmax, double outmin, double outmax){
		cvScale(this, this, (outmax-outmin)/(inmax-inmin), outmin-inmin*(outmax-outmin)/(inmax-inmin)); 
	}
	void convert_scale(double outmin, double outmax){
		double inmin,inmax;
		cvMinMaxLoc(this, &inmin, &inmax);
		convert_scale(inmin,inmax,outmin,outmax);
	}

	/** normalize the image data to between 255.0 and 0 */
	void normalize(){
		CvPoint minp,maxp;
		double min,max;
		cvMinMaxLoc(this, &min, &max, &minp, &maxp);
		if(max-min > 0){
			cvScale(this, this, 255.0/(max-min), 255*(-min)/(max-min));
		}
		else{
			cvScale(this, this, 1, -min);
		}
	}  
	
	const Image<T,nch> & operator= (const CvMat &mat){
		IplImage header;
		cvGetImage(&mat, &header);
		return (*this)=header;
	}

	/** copy operator for IplImage struct -- will handle naive format conversion (i.e. no scaling is done, only casting */
	const Image<T,nch> & operator= (const IplImage &im){
		assert(im.imageData!=NULL);
		if(im.imageData==NULL) return (*this);
		switch(im.depth){
			case IPL_DEPTH_8U:
				switch(im.nChannels){
					case 1:
						return this->convert(*(const Image<unsigned char, 1> *)(&im));
					case 2:
						return this->convert(*(const Image<unsigned char, 2> *)(&im));
					case 3:
						return this->convert(*(const Image<unsigned char, 3> *)(&im));
				}
				break;
			case IPL_DEPTH_16S:
				switch(im.nChannels){
					case 1:
						return this->convert(*(const Image<short, 1> *)(&im));
					case 2:
						return this->convert(*(const Image<short, 2> *)(&im));
					case 3:
						return this->convert(*(const Image<short, 3> *)(&im));
				}
				break;
			case IPL_DEPTH_32S:
				switch(im.nChannels){
					case 1:
						return this->convert(*(const Image<int, 1> *)(&im));
					case 2:
						return this->convert(*(const Image<int, 2> *)(&im));
					case 3:
						return this->convert(*(const Image<int, 3> *)(&im));
				}
				break;
			case IPL_DEPTH_32F:
				switch(im.nChannels){
					case 1:
						return this->convert(*(const Image<float, 1> *)(&im));
					case 2:
						return this->convert(*(const Image<float, 2> *)(&im));
					case 3:
						return this->convert(*(const Image<float, 3> *)(&im));
				}
				break;
			case IPL_DEPTH_64F:
				switch(im.nChannels){
					case 1:
						return this->convert(*(const Image<double, 1> *)(&im));
					case 2:
						return this->convert(*(const Image<double, 2> *)(&im));
					case 3:
						return this->convert(*(const Image<double, 3> *)(&im));
				}
				break;
		}
		return (*this);
	}

	/** return a pointer to the data at row i */
	T* operator[](const int &i) {
		return (T*) (this->imageData+(i*this->widthStep));
	}

	/** return a constant pointer to the data at row i */
	const T* operator[](const int &i) const {
		return (T*) (this->imageData+(i*this->widthStep));
	}
	
	/** return a constant pointer to the data at row i */
	const T* rowptr(const int &i) const {
		return (T*) (this->imageData+(i*this->widthStep));
	}

	/** return a pointer to the data at row i */
	T* rowptr(const int &i) {
		return (T*) (this->imageData+(i*this->widthStep));
	}

	/** allocate or reallocate (if neccessary) the image buffer.  
	    No memory allocation will take place if w==image.width and h==image.height
		@param w width of the reshaped image
		@param h height of the reshaped image
	*/
	bool reshape(int w, int h){
		std::cerr<<"cv::Image::reshape will go away soon"<<std::endl;
		return realloc( w, h);
	}
    bool realloc(const int & w, const int & h){
		if(w<=0 || h<=0) return false;
        if(this->imageData){
            if(this->width==w && 
               this->height==h){ 
                return false;
            }
			//printf("release(realloc): %p %d %d\n", this->imageData, this->width, this->height);
            cvReleaseData(this);
			this->imageData = NULL;
        }
		ImageHeader<T,nch>::_setParameters();
        cvInitImageHeader(this, cvSize(w,h), this->depth, this->nChannels);
        cvCreateData(this);
		//printf("create(realloc): %p %d %d\n", this->imageData, this->width, this->height);
		return true;
    }

	/** allocate or reallocate (if neccessary) the image buffer.
	  * @param im pointer to IplImage from which the width and height are taken
	  */
    bool reshape (const IplImage * im){
		std::cerr<<"cv::Image::reshape will go away soon"<<std::endl;
		return this->realloc(im);
	}
    bool realloc(const IplImage * im){
		if(im->roi){
			return this->realloc(im->roi->width, im->roi->height);
		}
		return this->realloc(im->width, im->height);
	}
	
	/** allocate or reallocate (if neccessary) the image buffer.
	  * @param im IplImage from which the width and height are taken
	  */
    bool reshape (const IplImage & im){
		std::cerr<<"cv::Image::reshape will go away soon"<<std::endl;
        return this->realloc(im);
    }
    bool realloc(const IplImage & im){
		if(im.roi){
			return this->realloc(im.roi->width, im.roi->height);
		}
        return this->realloc(im.width, im.height);
    }
    
	/** individual pixel access.  It is preferable to grab a pointer the whole row 
	    and iterate over that */
    inline const T & operator () (const int & x, const int & y, const int &ch=0) const {
        return *(T *)(this->imageData + y*this->widthStep + (x*this->nChannels + ch)*sizeof(T));    
        //return ((T *) this->imageData)[(y*this->width + x)*this->nChannels + ch];    
    }

	/** individual pixel access.  It is preferable to grab a pointer the whole row 
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
	const int & getDepth() const { return this->depth; }
	
	/** return the number of channels in the image */
	const int & getNumChannels() const { return this->nChannels; }

    /** use bilinear interpolation to calculate the pixel value at the given real number coordinates */
	T getSubPix(float fx, float fy) const{
		int x= cvFloor(fx), y=cvFloor(fy);
		int x2=x+1, y2=y+1;
		if(x2>=this->width) x2=x;
		if(y2>=this->height) y2=y;
		
		fx-=x;
		fy-=y;
		T ul,ur,ll,lr;
		ul = (*this)(x,y);
		ur = (*this)(x2,y);
		ll = (*this)(x,y2);
		lr = (*this)(x2,y2); 
		T ret =  (ul*(1-fx) + ur*fx) * (1-fy)   +
			     (ll*(1-fx) + lr*fx) * (fy);
		return ret;
    }
	
    /** */
	const Image<T,nch> & operator=(const T & val){
        CvScalar sc;
        for(int i=0; i<nch; i++){
            sc.val[i] = val;
        }
        cvSet(this, sc);
        return (*this);
    }
    
	/** subtract the given image from the original */
    const Image<T,nch> & operator -= (const Image<T,nch> & im){
        cvSub(this, &im, this);
        return (*this);
    }
	/** add the given image to the original */
    const Image<T,nch> & operator += (const Image<T,nch> & im){
        cvAdd(this, &im, this);
        return (*this);
    }
    
	/** multiply the elements of the original image by those of the argument */
    const Image<T,nch> & operator *= (const Image<T,nch> & im){
        cvMul(this, &im, this);
        return (*this);
    }
    
	/** return a new Image whose elements are the product of the elements of the arguments */
    Image<T,nch> operator * (const Image<T,nch> & im){
        Image<T,nch> out = (*this);
        return out*=im;
    }
    
    /*Image<T,nch> operator + (const Image<T,nch> & im){
        Image<T,nch> out = (*this);
        return out+=im;
    }
    
    Image<T,nch> operator - (const Image<T,nch> & im){
        Image<T,nch> out = (*this);
        return out-=im;
    }*/
    /*const Image<T,nch> & operator -= (const T val[nch]){
        CvScalar sc;
        for(int i=0; i<nch; i++){
            sc.val[i] = val[i];
        }
        cvSubS(this, sc, this);
        return (*this);
    }
    const Image<T,nch> & operator += (const T val[nch]){
        CvScalar sc;
        for(int i=0; i<nch; i++){
            sc.val[i] = val[i];
        }
        cvAddS(this, sc, this);
        return (*this);
    }
    
    const Image<T,nch> & operator *= (const T val[nch]){
        Image<T,nch> tim = (*this);
        tim = val;
        (*this)*=tim;
        return (*this);
    }
    
    Image<T,nch> operator * (const T val[nch]){
        Image<T,nch> out = (*this);
        return out*=val;
    }
    
    Image<T,nch> operator + (const T val [nch]){
        Image<T,nch> out = (*this);
        return out+=val;
    }
    
    Image<T,nch> operator - (const T val [nch] ){
        Image<T,nch> out = (*this);
        return out-=val;
    }*/
	/** swap the internal memory buffer of calling object and argument */
	void swap(Image<T,nch> & im){
		char * tmp;
		CV_SWAP(this->imageData, im.imageData, tmp);
		CV_SWAP(this->imageDataOrigin, im.imageDataOrigin, tmp);
		int itmp;
		CV_SWAP(this->width, im.width, itmp);
		CV_SWAP(this->height, im.height, itmp);
	}
	
	/** safe casting ensures that we don't try to cast an IplImage whose format does
	 *  not match this data type */
	static cv::Image<T, nch> * safe_cast( CvArr * arr ){

		int depth=cv::ImageDepth<T>::depth();
		IplImage * im;

		if(!CV_IS_IMAGE(arr)){
			fprintf(stderr, "Warning: Cannot cast non-IplImage to cv::Image\n");
			assert(0);
			return NULL;
		}

		im = (IplImage *) arr;
		if( depth != im->depth ){
			fprintf(stderr, "Warning: Trying to cast IplImage with depth=%d as %d\n", im->depth, depth);
			return NULL;
		}
		if( nch!=im->nChannels ){
			fprintf( stderr, "Warning: Trying to cast IplImage with nChannels=%d as cv::Image<T,%d>\n", im->nChannels, nch);
			return NULL;
		}
		return (cv::Image<T,nch> *) im;
	}
	static const cv::Image<T, nch> * safe_cast( const CvArr * arr ) {

        int depth=cv::ImageDepth<T>::depth();
        IplImage * im;

        if(!CV_IS_IMAGE(arr)){
            fprintf(stderr, "Warning: Cannot cast non-IplImage to cv::Image\n");
            assert(0);
            return NULL;
        }

        im = (IplImage *) arr;
        if( depth != im->depth ){
            fprintf(stderr, "Warning: Trying to cast IplImage with depth=%d as %d\n", im->depth, depth);
            return NULL;
        }
        if( nch!=im->nChannels ){
            fprintf( stderr, "Warning: Trying to cast IplImage with nChannels=%d as cv::Image<T,%d>\n", im->nChannels, nch);
            return NULL;
        }
        return (cv::Image<T,nch> *) im;
    }

		
	/** load the image data from the given filename 
	 *  this function supports gif,png,bmp,ppm,jpg,ipl formats */
	bool open(const char * fname){
		int len = strlen(fname);
		if(strcmp(fname+(len-3), "ipl")==0) return ipl_open(fname);
		if(strcmp(fname+(len-3), "ppm")==0) return pbm_open(fname);

		return cv_open(fname);
	}
	bool cv_open(const char * fname){

		IplImage * im = cvLoadImage(fname, (this->nChannels>1 ? 1 : 0));
		IplImage * me = (IplImage *) this;
		if(!im) return false;

		if(this->nChannels==im->nChannels && this->depth == im->depth){
			if(me->imageData || me->imageDataOrigin){
				//printf("release(open): %p %d %d\n", me->imageDataOrigin, me->width, me->height);
				cvReleaseData(me);
			}
			//memcpy(me, im, sizeof(IplImage));
			this->setData((T *) im->imageData, im->width, im->height);
			//printf("create(open): %p %d %d\n", me->imageDataOrigin, me->width, me->height);
			im->imageData=NULL;
			im->imageDataOrigin=NULL;
		}
		else{
			(*this) = (*im);
		}
		//printf("release(open): %p %d %d\n", im->imageDataOrigin, im->width, im->height);
		cvReleaseImage(&im);
		return true;
	}

	/** load the image data in ipl format from the given filename.
	 *  the ipl format is uncompressed and is provided mostly to provide
	 *  some functionality for saving image of type short,int,float or double 
     */
	bool ipl_open(const char * fname){
		// the binary flag is neccessary on windows!
		std::ifstream is(fname, std::ifstream::binary | std::ifstream::in);
		if(!is.is_open()) return false;

		bool ret = ipl_open(is);
		is.close();

		return ret;
	}
	bool pbm_open(const char * fname){
		char magic[2];
		FILE * f = fopen(fname, "r");
		int w,h,max;
		if(!f) return false;
		fread(magic, 1, 2, f);
		fscanf(f, "%d", &w);
		fscanf(f, "%d", &h);
		fscanf(f, "%d", &max);

		if(max > 255 && magic[1]=='6') { 
			fclose(f);
			return pbm_open_p6_16(fname, *this);	
		}
		else{
			fclose(f);
			return cv_open(fname);
		}

		return false;
	}
	template <typename data_t, int ch>
	static bool pbm_open_p6_16(const char * fname, cv::Image<data_t, ch> &im){
		cv::Image<unsigned short, 3> tmp;
		bool ret = pbm_open_p6_16(fname, tmp);
		im = tmp;
		return ret;
	}
	static bool pbm_open_p6_16(const char * fname, cv::Image<unsigned short, 3> &im){
		char magic[2];
		FILE * f = fopen(fname, "r");
		int w,h,max;
		if(!f) return false;
		fread(magic, 1, 2, f);
		fscanf(f, "%d", &w);
		fscanf(f, "%d", &h);
		fscanf(f, "%d", &max);

		if(max > 255 && magic[0]=='P' && magic[1]=='6') { 
			im.realloc(w, h);
			for(int i=0; i<h; i++){
				fread(im.imageData + im.widthStep*i, 3*sizeof(unsigned short), w, f);
			}
		}
		else{
			fclose(f);
			return im.cv_open(fname);
		}
		return true;
	}

	/** load the image data from the given stream.
	 *  @param is an input stream that has already been opened and is ready to be read
	*/
	bool ipl_open(std::istream & is){
		if(is.eof()) return false;
		
		IplImage header;
		unsigned long length;
		unsigned long offset;
		offset = is.tellg();
		//std::cout<<is.tellg()<<std::endl;
		is.seekg (0, std::ios::end);
		//std::cout<<is.tellg()<<std::endl;
		length = is.tellg();
		length-= offset;
		is.seekg (offset);
		//std::cout<<is.tellg()<<std::endl;	
		if(length<sizeof(IplImage)){
			std::cerr<<"Error reading IPL Header -- expected size of"<<sizeof(IplImage)<<" got "<<length<<std::endl;
			return false;
		}
		
		// read header
		is.read((char *) &header, sizeof(IplImage));
		
		//std::cout<<is.tellg()<<std::endl;	

		if(length<sizeof(IplImage)+header.imageSize){
			std::cerr<<"Erro reading Image data -- not large enough -- expected "<<sizeof(IplImage)+header.imageSize<<" got "<<length<<std::endl;
			return false;
		}
			
		// blank out pointer fields
		header.imageDataOrigin=0;
		header.imageId = 0;
		header.maskROI=0;
		header.roi=0;
		header.imageData=0;
		assert(header.imageSize!=0);
		
		if(header.depth == this->depth &&
		   header.nChannels == this->nChannels){

			// allocate memory if neccessary
			/*if(this->imageSize!=header.imageSize){
				if(this->imageData != NULL){
					IplImage * im = this;
					cvReleaseData(&im);
				}
				// simply copy header
				memcpy(this, &header, sizeof(IplImage));
				cvCreateData(this);
			}*/
			this->realloc(header);

			// read pixels
			is.read(this->imageData, this->imageSize);

		}
		else{
			cvCreateData(&header);
			//printf("create(ipl_open): %p %d %d\n", header.imageData, header.width, header.height);
			is.read(header.imageData, header.imageSize);
			
			// convert to my format
			(*this) = header;

			//printf("release(ipl_open): %p %d %d\n", this->imageData, this->width, this->height);
			cvReleaseData(&header);
		}
		//std::cout<<(*this)<<std::endl;
		return true;
	}
/*	friend std::ostream& operator << (std::ostream & os, const Image<T,nch> & im){
		//os<<"["<<im.width<<"x"<<im.height<<" - "<<nch<<" colors - "<<sizeof(T)<<" depth ]";

		CvMat m;
		cvGetMat(&im, &m);
		//cvReshape(&im,&m, nch, im.width*im.height);
		os<<*(Matrix<T> *) &m;
		return os;
	}
*/
	bool save(const char * fname) const{
		int len = strlen(fname);
		if(strcmp(fname+(len-3), "ipl")==0) return ipl_save(fname);
		if(strcmp(fname+(len-3), "pgm")==0) return pgm_save(fname);
		cvSaveImage(fname, this);
		return true;
	}
	bool ipl_save(const char *fname) const{
		// the binary flag is neccessary on windows!
		std::ofstream os(fname, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
		if(!os.is_open()) return false;
		bool ret = this->ipl_save(os);
		os.close();
		return ret;
	}	
	bool ipl_save(std::ostream & os) const{
		//std::cout<<os.tellp()<<std::endl;
		os.write((char *) this, sizeof(IplImage));
		//std::cout<<os.tellp()<<std::endl;
		os.write(this->imageData, this->imageSize);
		//std::cout<<os.tellp()<<std::endl;
		return true;
	}
	bool pgm_save(const char * fname) const{
		FILE * f = NULL;
		bool ret;
		f = fopen(fname, "w");
		if(f==NULL) return false;
		if((unsigned int) this->depth==IPL_DEPTH_16S || (unsigned int)this->depth==IPL_DEPTH_16U){
			ret = pgm_save_16(f);
		}
		else{
			ret = pgm_ascii_save(f);
		}
		fclose(f);
		return ret;
	}
	bool pgm_save_16(FILE * f) const{
		int nwrite;
		fprintf(f, "P5 %d %d 65535\n", this->width, this->height);
		if((nwrite = fwrite(this->imageData, 1, this->imageSize, f)) != this->imageSize){
			fprintf(stderr, "ERROR writing pgm, only %d of %d bytes written\n", nwrite, this->imageSize);
			return false;
		}
		return true;
	}
	bool pgm_ascii_save(FILE * f) const{
		fprintf(f, "P2\n");
		fprintf(f, "#Creator: cv::Image class\n");
		fprintf(f, "%d %d\n", this->width, this->height);
		fprintf(f, "256\n");
		
		double min,max,scale;
		cvMinMaxLoc(this, &min, &max);
		scale = 255.0/(max-min);
		for(int j=0; j<this->height; j++){
			const T * row = this->rowptr(j);
			for(int i=0; i<this->width; i++){
				fprintf(f, "%d ", (int)((row[i*nch]-min)*scale) );
				if(i%9==0) fprintf(f, "\n");
			}
			fprintf(f, "\n");
		}
		return true;
	}
	
	template <typename S, int nch2> 
	Image<T,nch> & eq_operator_same (const Image<S,nch2> &im){
		this->realloc(im);
		if(im.imageData==NULL || this->imageData==NULL) return (*this);
		cvCopy(&im, this);
		return (*this);
	}
	
	/** equals operator */
	Image<T,nch> & operator= (const Image<T,nch> & im){
		if(this->imageData==NULL && im.imageData==NULL) return (*this);
		assert(im.imageData!=NULL);
		this->realloc(im);
		cvCopy(&im, this);
		return (*this);
	}

	/** convert image of possibly type and color channels */
	template <typename S, int nch2> 
	Image<T,nch> & convert (const Image<S,nch2> & im){
		assert(im.imageData!=NULL);
		//if(im.imageData==NULL) return (*this) || this->imageData==NULL) return (*this);

		if(im.depth==this->depth){
			if(im.nChannels==this->nChannels){
				// this cast is to make the compiler happy 
				return this->eq_operator_same( im );
			}
			else{
				return this->eq_operator_nchannels( im );
			}
		}
		else if(im.nChannels==this->nChannels){
			return this->eq_operator_depth( im );
		}
		// need to convert twice here, so make sure to do less work
		// by doing type conversion on low color image
		else if(im.nChannels > this->nChannels){
			// convert color first
			Image<S,nch> dest;
			dest.eq_operator_nchannels(im);
			// now convert type
			return this->eq_operator_depth( dest );
		} 
		else{
			// convert type first
			Image<T,nch2> dest;
			dest.eq_operator_depth( im );
			// now convert color 
			return this->eq_operator_nchannels( dest );
		} 
	}

	/** convert image of one type to another -- color channels must be the same */
	template <typename S, int nch2> 
	Image<T,nch> & eq_operator_depth (const Image<S,nch2> & im){
		CvRect roi = cvGetImageROI(&im);
		this->realloc(im);
		assert(im.imageData!=NULL);
		assert(this->imageData!=NULL);
		for(int j=0; j<roi.height; j++){
			const S * s = im.rowptr(j+roi.y);
			T * t = this->rowptr(j);
			//std::cout<<"t="<<(unsigned long) t<<std::endl;
			//std::cout<<"s="<<(unsigned long) s<<std::endl;
			//ImageIterator<S> iit = im.row(j);
			//ImageIterator<T> tit = this->row(j);
			for(int i=0; i<roi.width*nch; i++){
				//tit[i] = (T) iit[i];
				//cvSet2D(this, j,i, cvGet2D(&im, j, i));
				t[i] = (T) (s[roi.x*3+i]);
				//std::cout<<"t[i]="<<(unsigned long) (&(t[i]))<<std::endl;
				//std::cout<<"s[i]="<<(unsigned long) (&(s[i]))<<std::endl;
			}
		}
		return (*this);
	}

	/** convert color to b/w and vice versa -- just knows about RGB for now */
	template <typename S, int nch2> 
	Image<T,nch> & eq_operator_nchannels (const Image<S,nch2> & im){
		this->realloc(im);
		if(im.imageData==NULL || this->imageData==NULL) return (*this);
		if(this->nChannels<im.nChannels && this->nChannels==1){
			if(this->depth==8){
				cvCvtColor(&im, this, CV_BGR2GRAY);
			}
			else{
				// I'm lazy, just use the first color channel as the conversion
				cvSplit(&im, this, NULL, NULL, NULL);
			}
		}
		else if(im.nChannels < this->nChannels && im.nChannels==1){
			for(int i=0; i<this->nChannels; i++){
				cvSetImageCOI(this, i+1);
				cvCopy(&im, this);
			}
			cvSetImageCOI(this,0);
		}
		return (*this);
	}
	
	Image<T,nch> dsample(int factor){
	}
	Image<T,nch> usample(int factor){
	
	}
	/** downsampes image by a factor of 2 */
	static Image<T,nch> & dscale2(const Image<T,nch> & a, Image<T,nch> & b) {
		CvRect roi = cvGetImageROI(&a);
		int w=roi.width,h=roi.height;
		b.realloc(w/2,h/2);
		//cvPyrDown(&a, &b);
		cvResize(&a, &b, CV_INTER_NN);
		return b;
	}

	static Image<T,nch> & uscale2(const Image<T,nch> & a, Image<T,nch> & b) {
		CvRect roi = cvGetImageROI(&a);
		int w=roi.width,h=roi.height;
		b.realloc(w*2,h*2);	
		cvPyrUp(&a,&b);
		return b;
	}
	
	// resize the image into the destination image 
	void resize(Image<T,nch> & dest, int method=CV_INTER_LINEAR) const{
		cvResize(this, &dest, method);
	}
	Image<T,nch> resize(int w, int h, int method=CV_INTER_LINEAR) const {
		Image<T,nch> im(w,h);
		cvResize(this, &im, method);
		return im;
	}
	
	Image<T,nch> scale(int w, int h, int method=CV_INTER_LINEAR) const {
		fprintf(stderr, "Scale has been renamed resize\n");
		assert(0);
		return (*this);
	}
	Image<T,nch> scale(float s, int method=CV_INTER_LINEAR) const{
		Image<T,nch> h((int)(this->width*s), (int)(this->height*s));
		for(int i=0; i<h.height; i++){
			for(int j=0; j<h.width; j++){
				h(j,i) = this->getSubPix(j/s, i/s);
			}
		}
		//cvResize(this, &h, method);
		return h;
	}
	void check(float v){
				assert(!isnan(v));
				assert(!isinf(v));
				assert(v<=255);
				assert(v>=0);
	}
	/*template <typename S, int n>
	void check(cv::Array<S,n> t){
		bool CANTm_CHECK_THIS_TYPE=false;
		assert(CANTm_CHECK_THIS_TYPE);
	}*/
	

	// downsample image
	/*Image<T,nch> dscale(int wdiv, int hdiv=-1, bool blur=true){
		Image<T,nch> h;
		if(hdiv<1) hdiv=wdiv;
		h.realloc(this->width/wdiv, this->height/hdiv);
		for(int i=0; i<h.height; i++){
			// iterator with step size of wdiv
			iterator it = this->row(i*hdiv, wdiv);
			iterator hit = h.row(i);
			for(size_t j=0; j<hit.length(); j++){
				hit[j] = it[j];	
			}
		}
		return h;
	}

	// upsample image
	Image<T,nch> uscale(int wmul, int hmul=-1){
		Image<T,nch> h;
		if(hmul<1) hmul = wmul;
		h.realloc(this->width*wmul, this->height*hmul);
		for(int i=0; i<h.height; i++){
			iterator it = this->row(i/hmul);
			iterator hit = h.row(i);
			for(int j=0; j<hit.length(); j++){
				hit[j] = it[j/wmul];
			}
		}
		return h;
	}*/
	Matrix<T, nch> asMatrix(){
		Matrix<T, nch> mat;
		cvGetMat(this, &mat);
		return mat;
	}
	~Image(){
		if(this->imageData) {
			//printf("release(realloc): %p %d %d\n", this->imageData, this->width, this->height);
			cvReleaseData(this);
		}
		if(this->roi){
			cvFree( (void **) (&(this->roi)) );
		}
        this->imageData = NULL;
	}

	// inefficient open mechanism
   	/*bool open(const char * fname){
		IplImage * im = cvLoadImage(fname);
		if(!im) return false;

		if(im->nChannels!=this->nChannels ||
		   im->depth!=this->depth){
		   	cvReleaseImage(&im);
			return false;
		}
		
		// free current image data
		if(this->imageData){
			cvReleaseData(this);
		}

		// copy over new data
		memcpy(this, im, sizeof(IplImage));

		// free old 
		cvReleaseImageHeader(&im);

		return true;
	}
	void show(const char * win){
		cvShowImage(win, this);
	}*/
	void show(const char * w) const{
		cvShowImage(w, this);
	}
	Image( const Image_OP<T,nch> & A){
		ImageHeader<T,nch>::_setParameters();
		this->realloc(A.getWidth(), A.getHeight());
		(*this) = A;
	}
	Image<T,nch> & operator = (const Image_OP<T,nch> & A){
		A.assign(*this);
		return (*this);
	}
	/*template <typename func_t>
	Image<T,nch> operator=(Array_BIN_OP<Image<T,nch>, func_t> op){
		op.func(*op.m_A, *op.m_B, (*this));
		return (*this);
	}*/
/*	Image<pixel_t, nch> & operator = (const Image_GEMM<pixel_t, nch> & A){
		A.assign(*this);
		return (*this);
	}
	Image<pixel_t, nch> & operator = (const Image_GEMM<pixel_t, nch> & A){
		A.assign(*this);
		return (*this);
	}
	Image<pixel_t, nch> & operator = (const Image_GEMM<pixel_t, nch> & A){
		A.assign(*this);
		return (*this);
	}*/
};


///////////////////////////// Image operations /////////////////////////////////////////////

template <typename pixel_t,int nch> struct Image_OP {
	virtual void assign(Image<pixel_t, nch> & X) const = 0;
	virtual int getWidth() const = 0;
	virtual int getHeight() const = 0;
};

template <typename pixel_t, int nch> struct Image_GEMM : Image_OP<pixel_t, nch> {
	typedef Image<pixel_t,nch> image_t;
	const image_t * m_A;
	const image_t * m_B;
	const image_t * m_C;
	double _alpha;
	double _beta;
	int _tABC;
	Image_GEMM(const image_t *A, const image_t *B, double alpha, const image_t * C, double beta, int tABC=0):
		m_A(A),
		m_B(B),
		m_C(C),
		_alpha(alpha),
		_beta(beta),
		_tABC(tABC)
	{
	}
	virtual void assign(Image<pixel_t,nch> & X) const{
		if(m_B!=NULL){
			cvGEMM(m_A, m_B, _alpha, m_C, _beta, &X, _tABC);
		}
	}
	virtual int getWidth() const {
		return m_A->roi ? m_A->roi.width : m_A->width;
	}
	virtual int getHeight() const {
		return m_A->roi ? m_A->roi.height : m_A->height;
	}
};

/// A+b
template <typename pixel_t, int nch> struct Image_ADDS : Image_OP<pixel_t, nch>{
	const Image<pixel_t,nch> * m_A;
	CvScalar _shift;
	Image_ADDS(const Image<pixel_t,nch> *A, CvScalar shift=0):
		        m_A(A),
				_shift(shift)
	{
	}
	virtual void assign(Image<pixel_t,nch> & X) const{
		cvAddS(m_A, _shift, &X);
	}
	virtual int getWidth() const {
		return m_A->width;
	}
	virtual int getHeight() const {
		return m_A->height;
	}

};

template <typename pixel_t, int nch> struct Image_ADD : Image_OP<pixel_t, nch>{
	const Image<pixel_t,nch> * m_A;
	const Image<pixel_t,nch> * m_B;
	double _scale1;
	double _scale2;
	double _shift;
	Image_ADD(const Image<pixel_t,nch> *A, const Image<pixel_t,nch> *B, double scale1, double scale2=1, double shift=0):
		m_A(A),
		m_B(B),
		_scale1(scale1),
		_scale2(scale2),
		_shift(shift)

	{
	}
	virtual void assign(Image<pixel_t,nch> & X) const{
		cvAddWeighted(m_A, _scale1, m_B, _scale2, _shift, &X);
	}
	virtual int getWidth() const {
		return m_A->width;
	}
	virtual int getHeight() const {
		return m_A->height;
	}
};

template <typename pixel_t, int nch> struct Image_SCALE : Image_OP<pixel_t, nch> {
	const Image<pixel_t,nch> * m_A;
	double _scale;
	Image_SCALE(const Image<pixel_t, nch> * A, double scale):
		m_A(A),
		_scale(scale){}
	virtual void assign(Image<pixel_t,nch> & X) const{
		cvScale(m_A, &X, _scale);
	}
	virtual int getWidth() const {
		return m_A->width;
	}
	virtual int getHeight() const {
		return m_A->height;
	}
};

template <typename pixel_t, int nch> struct Image_MUL : Image_GEMM<pixel_t, nch> {
	typedef Image<pixel_t,nch> image_t;
	Image_MUL(const image_t *A, const image_t *B):
		Image_GEMM<pixel_t,nch>(A,B,1.0,NULL,1.0)
	{
	}
};

// scale * A
template <typename pixel_t, int nch>
Image_SCALE<pixel_t,nch> operator * (const double & scale, const Image<pixel_t,nch> & A) {
	return Image_SCALE<pixel_t,nch>(&A, scale);
}
// A * scale
template <typename pixel_t, int nch>
Image_SCALE<pixel_t,nch> operator * (const Image<pixel_t,nch> & A, const double & scale) {
	return Image_SCALE<pixel_t,nch>(&A, scale);
}

// scale * (A* scale_n) == A*(scale*scale_n)
template <typename pixel_t, int nch>
Image_SCALE<pixel_t,nch> operator * (const double & scale, const Image_SCALE<pixel_t,nch> & A) {
	return Image_SCALE<pixel_t,nch>(A.m_A, A._scale*scale);
}

// (A* scale_n) * scale == A*(scale*scale_n)
template <typename pixel_t, int nch>
Image_SCALE<pixel_t,nch> operator * (const Image_SCALE<pixel_t,nch> & A, const double & scale) {
	return Image_SCALE<pixel_t,nch>(A.m_A, A._scale*scale);
}

// (A * scale) + B
template <typename pixel_t, int nch>
Image_ADD<pixel_t,nch> operator + (const Image_SCALE<pixel_t,nch> & A, const Image<pixel_t,nch> & B){
	return Image_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);
}

// A  + b
template <typename pixel_t, int nch>
Image_ADDS<pixel_t,nch> operator + (const Image<pixel_t,nch> & A, const CvScalar & shift){
	return Image_ADDS<pixel_t,nch>(&A, shift);
}

// A  - b
template <typename pixel_t, int nch>
Image_ADDS<pixel_t,nch> operator - (const Image<pixel_t,nch> & A, const CvScalar & shift){
	return Image_ADDS<pixel_t,nch>(&A, cvScalar(-shift.val[0], -shift.val[1], -shift.val[2], -shift.val[3]));
}

// A  + b
template <typename pixel_t, int nch>
Image_ADDS<pixel_t,nch> operator + (const Image<pixel_t,nch> & A, const double & shift){
	return Image_ADDS<pixel_t,nch>(&A, cvScalarAll(shift));
}

// A  - b
template <typename pixel_t, int nch>
Image_ADDS<pixel_t,nch> operator - (const Image<pixel_t,nch> & A, const double & shift){
	return Image_ADDS<pixel_t,nch>(&A, cvScalarAll(-shift));
}

// B + (A * scale)
template <typename pixel_t, int nch>
Image_ADD<pixel_t,nch> operator + (const Image<pixel_t,nch> & B, const Image_SCALE<pixel_t,nch> & A){
	return Image_ADD<pixel_t,nch>(A.m_A, &B, A._scale, 1);
}
// (A*scale) + (B * scale)
template <typename pixel_t, int nch>
Image_ADD<pixel_t,nch> operator + (const Image_SCALE<pixel_t,nch> & A, const Image_SCALE<pixel_t,nch> & B){
	return Image_ADD<pixel_t,nch>(A.m_A, B.m_A, A._scale, B._scale);
}

// (A+B) * scale
template <typename pixel_t, int nch>
Image_ADD<pixel_t,nch> operator * (const Image_ADD<pixel_t,nch> & A, const double & scale){
	return Image_ADD<pixel_t,nch>(A.m_A, A.m_B, A._scale1*scale, A._scale2*scale);
}

// A*B
template <typename pixel_t, int nch>
Image_MUL<pixel_t,nch> operator * (const Image<pixel_t,nch> & A, const Image<pixel_t,nch> &B) {
	return Image_MUL<pixel_t,nch>(&A, &B);
}

// A*B * scale
template <typename pixel_t, int nch>
Image_GEMM<pixel_t,nch> operator * (const Image<pixel_t,nch> & A, const Image_SCALE<pixel_t,nch> & B) {
	return Image_GEMM<pixel_t,nch>(A.m_A, B.m_A, B._scale, NULL, 0);
}

// A*B * scale
template <typename pixel_t, int nch>
Image_GEMM<pixel_t,nch> operator * (const Image_MUL<pixel_t,nch> & A, const double & scale) {
	return Image_GEMM<pixel_t,nch>(A.m_A, A.m_B, scale, NULL, 0);
}

// A+B
template <typename pixel_t,int nch>
Image_ADD<pixel_t, nch> operator + (const Image<pixel_t,nch> & A, const Image<pixel_t,nch> &B) {
	return Image_ADD<pixel_t,nch>(&A, &B, 1.0, 1.0);
}

// A-B
template <typename pixel_t, int nch>
Image_ADD<pixel_t, nch> operator - (const Image<pixel_t, nch> & A, const Image<pixel_t, nch> &B) {
	return Image_ADD<pixel_t, nch>(&A, &B, 1.0, -1.0);
}

// A*B+C
template <typename pixel_t, int nch>
Image_GEMM<pixel_t, nch> operator + (Image_MUL<pixel_t, nch> & M, const Image<pixel_t, nch> &C){
	return Image_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);
}

//C+A*B
template <typename pixel_t, int nch>
Image_GEMM<pixel_t, nch> operator + (const Image<pixel_t, nch> &C, Image_MUL<pixel_t, nch> & M){
	return Image_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, 1.0);
}

//C - A*B
template <typename pixel_t, int nch>
Image_GEMM<pixel_t, nch> operator - (const Image<pixel_t, nch> &C, Image_MUL<pixel_t, nch> & M){
	return Image_GEMM<pixel_t, nch>(M.m_A, M.m_B, -1.0, &C, 1.0);
}

//A*B - C
template <typename pixel_t, int nch>
Image_GEMM<pixel_t, nch> operator - (Image_MUL<pixel_t, nch> & M, const Image<pixel_t, nch> &C){
	return Image_GEMM<pixel_t, nch>(M.m_A, M.m_B, 1.0, &C, -1.0);
}

} // namespace cv

#if 0
int main(){
    cv::Image<int ,3> im1(39,39);
    cv::Image<int ,3> im(1,1);
    cv::Image<int ,1> redch(1,1);
    redch = im1.getChannel(0);
    /*for(int i=0; i<39; i++){
        for(int j=0; j<39; j++){
            im1(i,j,0) = i;
            im1(i,j,1) = j;
            im1(i,j,2) = i*j;
        }
    }
    */
    im = im1;
    int val[3] = {5, 4, 2};
    im = val;
    /*im-=val;
    im+=val;
    */im*=val;//val;
    /*for(int i=0; i<39; i++){
        for(int j=0; j<39; j++){
            printf("%d,%d, %d = %d\n", i, j, 0, im(i,j,0));
   //         assert(im(i,j,0)==cvGet2D(&im, j, i).val[0]);
   //         assert(im(i,j,0)==i);
            printf("%d,%d, %d = %d\n", i, j, 1, im(i,j,1));
   //         assert(im(i,j,1)==cvGet2D(&im, j, i).val[1]);
    //        assert(im(i,j,1)==j);
            printf("%d,%d, %d = %d\n", i, j, 2, im(i,j,2));
            //assert(im(i,j,2)==cvGet2D(&im, j, i).val[2]);
            //assert(im(i,j,2)==i*j);
        }
    }*/
    
}
#endif

#endif //__CV_IMAGE_HH
