#include <cv/image.h>
#include <cv/dog.h>
// extract individual ccd pixel responses from bayer ccd 
// ___________
// | q1 | q2 |
// | -  + -  |
// | q3 | q4 |
// -----------
template <typename data_t>
void cvBayerSplit(const cv::Image<data_t, 1> & im, 
		          cv::Image<data_t, 1> & q1,
		          cv::Image<data_t, 1> & q2,
		          cv::Image<data_t, 1> & q3,
		          cv::Image<data_t, 1> & q4){
	q1.reshape(im.width/2, im.height/2);
	q2.reshape(im.width/2, im.height/2);
	q3.reshape(im.width/2, im.height/2);
	q4.reshape(im.width/2, im.height/2);

	for(int i=0; i<im.height; i++){
		const data_t * p = im[i];
		data_t * pq1, *pq2;
		if(i%2==0){
			pq1=q1[i/2];
			pq2=q2[i/2];
		}
		else{
			pq1=q3[i/2];
			pq2=q4[i/2];
		}
		for(int j=0; j<im.width; j+=2){
			pq1[j/2] = p[j];
			pq2[j/2] = p[j+1];
		}
	}
}
int cvSegReport( int status, const char * func_name,
		               const char *err_msg, const char * file_name,
					                  int line, void *userdata ){
	fprintf(stderr, "ERROR in function %s at file %s:%d:\n%s\n", func_name, file_name, line, err_msg);
	    assert(0);
}

int main(int argc, char ** argv){
	cv::Image<float, 1> im, r, g, b, e;
	cvRedirectError(cvSegReport);
	assert(im.open(argv[1]));
	cvBayerSplit(im,  r,g,b,e);

	b.normalize();
	cv::Image<uchar> x = b;
	cv::Image<float> IR = b;
	cv::Image<float> IR2 = b;
	cv::Image<uchar> bin= b;
	cvNamedWindow("win", 0);
	x.show("win");
	cvWaitKey(-1);
	cvAdaptiveThreshold(&x, &bin, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 5, -110 );
	
	CvPoint min,max1,max2;
	double minv,maxv;
	cvMinMaxLoc(&bin, &minv, &maxv, &min, &max1);
	cvFloodFill(&bin, max1, CV_RGB(0,0,0));
	cvMinMaxLoc(&bin, &minv, &maxv,  &min, &max2);

	r.normalize();
	cv::Image<uchar,3> im3 = r;

	cvCircle(&im3, max1, 3, CV_RGB(0,255,0));
	cvCircle(&im3, max2, 3, CV_RGB(0,255,0));
	im3.show("win");
	cvWaitKey(-1);

	//r = r.scale(500,500);
	//b = b.scale(500,500);
	//g = g.scale(500,500);
	//e = e.scale(500,500);
	/*cvNamedWindow("IR",0);
	IR.imagesc("IR");
	cvWaitKey(-1);
	IR2.imagesc("IR");
	cvWaitKey(-1);
	b.imagesc("IR");
	cvWaitKey(-1);*/
	return 0;
}
