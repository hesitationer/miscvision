#include <highgui.h>
#include <cv/image.hh>

template <typename T>
void smooth(const cv::Image<int> & integral, int radius, cv::Image<T> & result){
	int w,h;
	w = result.width-radius;
	h = result.height-radius;
	cvZero(&result);
	for(int j=radius; j<h; j++){
		for(int i=radius; i<w; i++){
			result[j][i] = (integral[j+radius+1][i+radius+1] + integral[j-radius][i-radius] -
			                integral[j-radius][i+radius+1] - integral[j+radius+1][i-radius]) / 
							   ((2*radius+1)*(2*radius+1));
		}
	}
	int offset = w;
	int I,J;
	for(int j=radius; j<h; j++){
		for(int i=0; i<radius; i++){
			result[j][i] = (integral[j+radius+1][i+radius+1] + integral[j-radius][0] -
			                integral[j-radius][i+radius+1] - integral[j+radius+1][0]) / 
							   ((radius+i+1)*(2*radius+1));
			I = i+w;
			result[j][I] = (integral[j+radius+1][integral.width-1] + integral[j-radius][I-radius] -
					integral[j-radius][integral.width-1] - integral[j+radius+1][I-radius])/ 
				((2*radius - i)*(2*radius+1));
		}
	}
	
	for(int i=radius; i<w; i++){
		for(int j=0; j<radius; j++){
			result[j][i] = (integral[j+radius+1][i+radius+1] + integral[0][i-radius] -
			                integral[j+radius+1][i-radius] - integral[0][i+radius+1]) / 
							   ((radius+j+1)*(2*radius+1));
			J = j+h;
			result[J][i] = (integral[integral.height-1][i+radius+1] + integral[J-radius][i-radius] -
					integral[integral.height-1][i-radius] - integral[J-radius][i+radius+1])/ 
				((radius + 1 + (radius-j-1))*(2*radius+1));
		}
	}
	for(int j=0; j<radius; j++){
		for(int i=0; i<radius; i++){
			result[j][i] = (integral[j+radius+1][i+radius+1] + integral[0][0] -
			                integral[j+radius+1][0] - integral[0][i+radius+1]) / 
							   ((radius+i+1)*(radius+j+1));
			I = i+w;
			J = j+h;
			result[j][I] = (integral[j+radius+1][integral.width-1] + integral[0][I-radius] -
					integral[0][integral.width-1] - integral[j+radius+1][I-radius])/ 
				((2*radius - i)*(radius+j+1));
			result[J][i] = (integral[integral.height-1][i+radius+1] + integral[J-radius][0] -
					integral[integral.height-1][0] - integral[J-radius][i+radius+1])/ 
				((radius + 1 + (radius-j-1))*(radius+i+1));
			result[J][I] = (integral[integral.height-1][integral.width-1] + integral[J-radius][I-radius] -
					integral[integral.height-1][I-radius] - integral[J-radius][integral.width-1])/ 
				((radius + 1 + (radius-j-1))*(2* radius - i ));
				
		}
	}
}
int count_lt(cv::Image<float> &im, const float & v, int x, int y){
	int ct=0;
	for(int ii=-1; ii<=1; ii++){
		for(int jj=-1; jj<=1; jj++){
			if(im(x+ii,y+jj)<v) ct++;
		}
	}
	return ct;
}
int count_gt(cv::Image<float> &im, const float & v, int x, int y){
	int ct=0;
	for(int ii=-1; ii<=1; ii++){
		for(int jj=-1; jj<=1; jj++){
			if(im(x+ii,y+jj)>v) ct++;
		}
	}
	return ct;
}
void count_gt(cv::Image<float> & before, cv::Image<float> & current, cv::Image<float> & after, cv::Image<uchar> & count, cv::Image<uchar> & minim, cv::Image<uchar> & maxim){
	float v;
	float ct;
	cvZero(&minim);
	cvZero(&maxim);
	
	for(int j=1; j<current.height-1; j++){
		for(int i=1; i<current.width-1; i++){
			ct = 0;
			v = current(i,j);
			maxim(i,j) = count_gt(before, v, i, j) +
			             count_gt(current, v, i, j) +
						 count_gt(after, v, i, j);
			minim(i,j) = count_lt(before, v, i, j) +
			             count_lt(current, v, i, j) +
						 count_lt(after, v, i, j);
		}
	}
	CvPoint min,max;
	double minv,maxv;
	//cvMinMaxLoc(&count, &minv, &maxv, &min, &max);
	std::cout<<"min = "<<minv<<std::endl;
	std::cout<<"max = "<<maxv<<std::endl;
	cvThreshold(&maxim, &maxim, 25, 255, CV_THRESH_BINARY);
	cvThreshold(&minim, &minim, 25, 255, CV_THRESH_BINARY);
	cvOr(&maxim,&minim, &count);
}

int main(int argc, char** argv){
	cv::Image<uchar,3> im3,cp3;
	cv::Image<int> integral;
	cv::Image<float> sm1;
	cv::Image<float> sm2;
	cv::Image<float> sm3;
	cv::Image<uchar> im,count,min,max;
	cv::Image<float> f1,f2,d1,d2,d3;
	
	cvNamedWindow("smoothed", 0);	
	for(int i=1; i<argc; i++){
		if(!im3.open(argv[i])) return -1;
		im = im3;
		sm2 = im;
		im.show("smoothed");
		cvWaitKey(1);
		integral.reshape(im.width+1, im.height+1);
		sm1.reshape(sm2);
		min.reshape(im);
		max.reshape(im);
		count.reshape(im);
		d1 = sm2;
		d2 = sm2;
		d3 = sm2;
		cvIntegral(&im, &integral, NULL,NULL);
		float scale = 2.0;
		for(float j=1; j<200; j+=scale){
			scale*=2.0;
			std::cout<<j<<std::endl;
			smooth(integral, j, sm1);
			d1 = d2;
			d2 = d3;
			d3 = sm2 - sm1;
			sm2 = sm1;
			count_gt(d1,d2,d3,count,min,max);
			cp3 = im3;
			for(int j=0; j<count.height; j++){
				for(int i=0; i<count.width; i++){
					if(count(i,j)>0){
						cvCircle(&cp3, cvPoint(i,j), 2, CV_RGB(255,0,0));
					}
				}
			}
			cp3.show("smoothed");
			cvWaitKey(1);
		}
	}
		
}
