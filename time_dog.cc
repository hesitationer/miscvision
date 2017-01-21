#include <cv/gaussian.hh>
#include <cv/image.hh>
#include <highgui.h>

class gcoef {
public: 
	int len;
	Gaussian g;
	float * coefs;
	gcoef(float s, int l): len(l), g(l/2+1,s), coefs(new float[len]) 
	{
		for(int i=0; i<len; i++){
			coefs[i] = g.pdf(i);
		}
	}
};

int main(int argc, char ** argv){
	gcoef gc[4] = { gcoef(1.0, 21), gcoef(2.0, 21), gcoef(3.0, 21), gcoef(6.0, 21) };
	cv::Image<uchar> im;
	cv::Image<float> g[4];
	cv::Image<float> d[4];
	cv::Image<float> f[21];
	cv::Image<uchar> cmp1,cmp2;
	char * windowg[4]  = { "g1", "g2", "g3", "g4" };
	char * windowd[4]  = { "d1", "d2", "d3", "d4" };
	for(int i=0; i<4; i++){
		cvNamedWindow(windowg[i], 0);
		cvNamedWindow(windowd[i], 0);
	}
	for(int i=0; i<21; i++){
		im.open(argv[i+1]);
		f[i] = im;
		cmp1.reshape(im);
		cmp2.reshape(im);
	}
	for(int i=21; i<(argc-1); i++){
		im.open(argv[i+1]);
		f[i%21] = im;		
		
		for(int j=0; j<4; j++){
			g[j].reshape(f[0]);
			cvZero(&(g[j]));
			d[j].reshape(f[0]);
		}
		// calculate gaussians
		for(int k=0; k<4; k++){
			for(int j=0; j<21; j++){
				d[k] = g[k]+f[j]*gc[k].coefs[j];
				g[k] = d[k];
			}
			//g[k].convert(1.0,0).show(windowg[k]);
			cvWaitKey(1);
		}
		for(int k=1; k<4; k++){
			d[k-1] = g[k-1] - g[k];
			d[k-1].convert(.5, 128).show(windowd[k]);
			cvWaitKey(1);
		}
		cvCmp(&(d[1]), &(d[0]), &cmp1, CV_CMP_LT);
		cvCmp(&(d[1]), &(d[2]), &cmp2, CV_CMP_LT);
		cvAnd(&cmp1, &cmp2, &cmp1);
		cmp1.show(windowg[0]);
		
		cvCmp(&(d[1]), &(d[0]), &cmp1, CV_CMP_GT);
		//cvCmp(&(d[1]), &(d[2]), &cmp2, CV_CMP_GT);
		//cvAnd(&cmp1, &cmp2, &cmp1);
		cmp1.show(windowg[1]);
		cvWaitKey(1);
		
	}
}
