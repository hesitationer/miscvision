#ifndef __COLOR_MAP_HH
#define __COLOR_MAP_HH

#include <vector>
#include <ostream>
#include <cv/matrix.hh>
#include <cv/image.h>
#include <cv/array.h>

// maps indices to colors
class ColorMap : public std::vector< cv::Array<uchar, 3> > {
public:
	typedef cv::Array<uchar, 3> pixel_t;	
	typedef std::vector< pixel_t > base_t;	
	typedef enum {
		Grey,
		Hot,
		Ocean,
		Jet
	} ColorMapEnum;
	ColorMap(){
	}
	ColorMap(int size) : base_t(size){
	}
	ColorMap(ColorMapEnum type){
		(*this)(type);
	}
	CvScalar lookup(int i){
		i = MAX(MIN(i, this->size()-1), 0);
		cv::Array<uchar,3> c = (*this)[i];
		return cvScalar(c[2],c[1],c[0]);
	}
	void operator()(ColorMapEnum type){
		switch(type){
		case Ocean:
			this->OceanMap(255);
			break;
		case Jet:
			this->JetMap(255);
			break;
		case Hot:
			this->HotMap(255);
			break;
		case Grey:
		default:
			this->GreyMap(255);
			break;
		}
	}
	void GreyMap(int levels){
		if(levels>255) levels=255;
		this->resize(levels);
		for(int i=0; i<levels; i++){
			uchar c = (i*255/levels);
			(*this)[i]=pixel_t(c,c,c);
		}
	}
	void HotMap(int levels){
		double x;
		if(levels<0) levels=255;
		this->resize(levels);
		for(int i=0; i<levels; i++){
			x=(double)i/(levels-1);
			// B,G,R
			(*this)[i][0] = (uchar) (((x >= 4./5.) * (5.*x - 4.)) * 255);  
			(*this)[i][1] = (uchar) (((x >= 2./5. & x < 4./5.) * (5./2. * x - 1.) + (x >= 4./5.)) * 255);
			(*this)[i][2] = (uchar) (((x < 2./5.) * (5./2. * x) + (x >= 2./5.)) * 255);
		}	
	}
	void OceanMap(int levels){
		float cutin = (float)(levels/3);
		if(levels<0) levels=255;
		this->resize(levels);
		for(int i=0; i<levels; i++){
    		(*this)[i][0] = (uchar) (i*255/(levels-1));
			(*this)[i][1] = (uchar) (i<cutin ? 0 : (i-cutin)/(levels-cutin-1)*255);
			(*this)[i][2] = (uchar) (i<(2*cutin) ? 0 : (i-2*cutin)/(levels-2*cutin-1)*255);
		}
	}
	void JetMap(int levels){
		if(levels<0) levels=255;
		this->resize(levels);
		float x,r,g,b;
		for(int i=0; i<levels; i++){
			x = (float)i/(levels-1);
	    	r = (x >= 3./8. & x < 5./8.) * (4. * x - 3./2.)\
				+ (x >= 5./8. & x < 7./8.) + (x >= 7./8.) * (-4. * x + 9./2.);
			g = (x >= 1./8. & x < 3./8.) * (4. * x - 1./2.)\
				+ (x >= 3./8. & x < 5./8.) + (x >= 5./8. & x < 7./8.) * (-4. * x + 7./2.);
			b = (x < 1./8.) * (4. * x + 1./2.) + (x >= 1./8. & x < 3./8.)\
				+ (x >= 3./8. & x < 5./8.) * (-4. * x + 5./2.);
			(*this)[i][0] = (uchar) (b*255);
			(*this)[i][1] = (uchar) (g*255);
			(*this)[i][2] = (uchar) (r*255);
		}
	}
	void RedMap(int levels){
		if(levels<0) levels=255;
		this->resize(levels);
		float x,r,g,b;
		for(int i=0; i<levels/2; i++){
			r = ((float)i)/(levels/2);
			(*this)[i][0] = (uchar) 0; 
			(*this)[i][1] = (uchar) 0; 
			(*this)[i][2] = (uchar) (r*255);
		}
		for(int i=levels/2; i<levels; i++){
			g = b = ((float)(i-levels/2))/(levels/2);
			(*this)[i][0] = (uchar) (b*255);
			(*this)[i][1] = (uchar) (g*255);
			(*this)[i][2] = (uchar) 255;
		}
	}
	template <typename spixel_t>
	cv::Image<uchar,3> & remap(const cv::Image<spixel_t,1> & src, cv::Image<uchar,3> & dest, double min=-1, double max=-1) const{
		int size = this->size()-1;
		int idx;
		
		// find min/max of src if unspecified
		if(min==-1 && max==-1){
			cvMinMaxLoc(&src, &min, &max);
		}
		else if(min==-1){
			double t;
			cvMinMaxLoc(&src, &min, &t);
		}
		else if(max==-1){
			double t;
			cvMinMaxLoc(&src, &t, &max);
		}
		
		for(int j=0; j<src.height; j++){
			const spixel_t * srow = src[j];
			pixel_t * drow = (pixel_t *) dest[j];
			for(int i=0; i<src.width; i++){
				idx = MAX(0, MIN(size-1, (int)((size * (srow[i] - min))/(max-min))));
				drow[i] = (*this)[idx];
			}
		}
		return dest;
	}
	template <typename spixel_t>
	void remap(const Matrix<spixel_t> & src, cv::Image<uchar, 3> & dest, double min=-1, double max=-1) const{
		int size = this->size()-1;
		int idx;
		
		// find min/max of src if unspecified
		if(min==-1 && max==-1){
			cvMinMaxLoc(&src, &min, &max);
		}
		else if(min==-1){
			double t;
			cvMinMaxLoc(&src, &min, &t);
		}
		else if(max==-1){
			double t;
			cvMinMaxLoc(&src, &t, &max);
		}
		
		for(int j=0; j<src.height; j++){
			const spixel_t * srow = src[j];
			pixel_t * drow = (pixel_t *) dest[j];
			for(int i=0; i<src.width; i++){
				idx = MAX(0, MIN(size-1, (int)((size * (srow[i] - min))/(max-min))));
				drow[i] = (*this)[idx];
			}
		}
		return dest;
	}
	friend std::ostream & operator<<(std::ostream & out, const ColorMap & c){
		for(size_t i=0; i<c.size(); i++){
			out<<c[i]<<std::endl;
		}
		return out;
	}
};

template< typename T, int nch>
void imagesc(const char * window, const cv::Image<T,nch> & im, const ColorMap & m, int min=-1, int max=-1){
	cv::Image<uchar, 3> im3(im.width, im.height);
	m.remap(im, im3, min, max).show(window);
}


#endif
