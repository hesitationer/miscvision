#ifndef __DEPTH_IMAGE_HH
#define __DEPTH_IMAGE_HH

#include <cv/image.h>
#include <cv/vector3.h>

namespace cv {

template
<typename i_type=uchar, typename z_type=float, int nchannels=3>
class DepthImage {
protected:
	double _cx;
	double _cy;
	double _centerX;
	double _centerY;
public:
	// the image are public for easy modification
	cv::Image<i_type, nchannels> intensity;
	cv::Image<z_type> depth;
	
	DepthImage() {}
	DepthImage(double cx, double cy, double cenX, double cenY) :
		_cx(cx),_cy(cy),_centerX(cenX), _centerY(cenY){
	}
	fVec3 to3(const double &u, const double &v, const double &z){
		return fVec3((u-_centerX)*_cx*z, (v-_centerY)*_cy*z, z);
	}
	fVec3 to3(const float &u, const float &v, const float &z){
		return fVec3((u-_centerX)*_cx*z, (v-_centerY)*_cy*z, z);
	}
	fVec3 to3(const double &u, const double &v) {
		fVec3 X;
		X[2] = depth[(int)v][(int)u];
		X[0] = (u-_centerX)*_cx*X[2];
		X[1] = (v-_centerY)*_cy*X[2];
		return X;
	}
	fVec3 to3(const fVec2 & x){
		return to3(x[0],x[1]);
	}
	fVec2 to2(const double &x, const double &y, const double &z){
		fVec2 X;
		X[0] = x/(z*_cx) + _centerX;
		X[1] = y/(z*_cy) + _centerY;
		return X;
	}
	fVec2 to2(const fVec3 & X){
		return to2(X[0], X[1], X[2]);
	}
	const double & getCx() { return _cx; }
	const double & getCy() { return _cy; }
	const double & getZCenterU() { return _centerX; }
	const double & getZCenterV() { return _centerY; }
	void setCx(const double & cx) { _cx=cx; }
	void setCy(const double & cy) { _cy=cy; }
	void setZCenterU(const double & cenX) { _centerX = cenX; }
	void setZCenterV(const double & cenY) { _centerY = cenY; }
	bool isValid(const int & u, const int &v){
		return u>=0 && v>=0 && u<depth.width && v<depth.height && depth[v][u]!=0;
	}
};
	
} // namespace cv

#endif // __DEPTH_IMAGE
