#ifndef __FILTERED_IMAGE_HH
#define __FILTERED_IMAGE_HH

#include "cv/image.h"

/**
 * Base class for image filters.  The default usage is to apply a filter to an 
 * existing image.  This filter then contains an image data that can be used 
 * directly or copied. 
 * @author Roman Stanchak
 */
class FilteredImage : public RImage {
public:
    ///process input image using filter
	void filter(IplImage * im, bool in_place=true){
        this->myFilter(im, in_place);
    }
    ///process input image using filter
	void filter(RImage & r, bool in_place=true){
        this->myFilter(r.getIplImage(), in_place);
    }
protected: 
    /// implementation method for filter
    /** 
    * method containing the implementation for the filter.  By default the filter
    * will copy the input image.  When subclassing, overriding this method is 
    * preferable
    * @param im input image.
    * @param in_place specifies whether to operate on im or allocate new data
    */
    virtual void myFilter(IplImage *im, bool in_place=true) = 0;
	/*{
		if(in_place){
			this->setIplImage(im);
		}	
		else{
			RImage::reshape(& m_image, im);
			cvCopy(im, m_image);
		}
	}*/
    
    
/*    virtual RImage * applyFilter(FilteredImage * f, RImage * out){
        bool in_place = (out==NULL);
        f->filter(m_image, in_place);
        if(!in_place){
            out->reshape(f);
            cvCopy(f->getIplImage(), out->getIplImage());
            return out;
        }
        return this;

    }
*/
};

#endif //__FILTERED_IMAGE_HH
