#include <cv/image.h>
#include <highgui.h>
#include <iostream>

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
	cvMinMaxLoc(&maxim, &minv, &maxv, &min, &max);
	std::cout<<"min = "<<minv<<std::endl;
	std::cout<<"max = "<<maxv<<std::endl;
	cvThreshold(&maxim, &maxim, 18, 255, CV_THRESH_BINARY);
	cvThreshold(&minim, &minim, 18, 255, CV_THRESH_BINARY);
	cvOr(&maxim,&minim, &count);
}


int main(int argc, char ** argv){
	cv::Image<uchar> im,count,min,max;
	cv::Image<uchar, 3 > im3, cp3;
	im3.open(argv[1]);
	cv::Image<float> f1,f2,d1,d2,d3;
	int gwidth = 11;
	double sigma = (3*.5-1.0)*.3 + .8;

	//im3 = im3.scale(100,100);
	im = im3;
	f1 = im;
	count = im;
	min = im;
	max = im;
	d3 = f1;
	d2 = d3;
	d1 = d2;
	
	cvNamedWindow("w", 0);
	cvNamedWindow("d1", 0);
	cvNamedWindow("d2", 0);
	cvNamedWindow("d3", 0);
	for(int i=0; i<1000; i++){
		d3 = d2;
		d2 = d1;
		f2 = f1;
		
		cvSmooth(&f1, &f1, CV_GAUSSIAN, gwidth, gwidth);
		d1 = f1 - f2;
		cvAbs(&d1, &d3);
		/*count_gt(d1,d2,d3,count,min,max);
		cp3=im3;
		for(int j=0; j<count.height; j++){
			for(int k=0; k<count.width; k++){
				if(count(k,j)>0){
					cvCircle(&cp3, cvPoint(k,j), sigma*(i+1), CV_RGB(255,0,0), 1, CV_AA);
				}
			}
		}*/
		count.show("w");
		d1.imagesc("d1");
		d2.imagesc("d2");
		d3.imagesc("d3");
		std::cout<<i<<std::endl;
		if(cvWaitKey(-1)=='q') break;
	}	
}
