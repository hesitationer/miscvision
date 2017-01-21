#include <fixed_priority_queue.hpp>
#include <cvarr_adapter.hpp>

template<typename T1, int nch1, typename T2, int nch2>
int cvCopyN( cv::Image<T1, nch1> & im1, cv::Image<T2, nch2> & im2){
	printf("cvCopyN called\n");
	return 0;
}
template<typename T1, int nch1>
int cvCopyN( cv::Image<T1, nch1> & im1, cv::Image<double, 1> & im2){
	printf("cvCopyN not supported for double\n");
	return 0;
}

CVARR_TO_IMAGE_FUNC2( cvCopyN )

struct cvTopNAdapter {
	int num;
	float min; 
	float max;
	bool invert;

	cvTopNAdapter() :
		num(50),
		invert(false)
	{}

	template<typename T, int nch>
	int operator() ( cv::Image<T,nch> & im ){
		printf("cvTopN not implemented for channels > 1\n");
		return -1;
	}

	template<typename T>
	int operator() ( cv::Image<T,1> & im ){
		fixed_priority_queue< CvPoint > queue(num);
		float scale = (invert ? -1.0 : 1.0);
		CvRect roi = im.getImageROI();
		roi.width+=roi.x;
		roi.height+=roi.y;
		for(int j=roi.y; j<roi.height; j++){
			for(int i=roi.x; i<roi.width; i++){
				//printf("%f  ?=  ", im[j][i]);
				//print_queue(queue);
				queue.push(cvPoint(i,j), scale*im[j][i]);
			}
		}
		//printf("----------\n");
		min = queue.top().weight;
		queue.pop();
		while(!queue.empty()){
			max = queue.top().weight;
			//printf("%f\n", min);
			queue.pop();
		}
		return 0;
	}
};

void cvTopN( CvArr * arr, int num, float * min, float * max ){
	cvTopNAdapter base_func;
	base_func.num = num;
	IplImageAdapter1< cvTopNAdapter > func( base_func );
	func( arr );
	*min = base_func.min;
	*max = base_func.max;
}

template <typename data_t, int nch> 
int cvPrintit( cv::Image<data_t, nch>  & im ){
	printf("cv::Image<%ld,%d>\n", sizeof(data_t), nch);
	return 0;
}
int cvPrintit( cv::Image<double, 3>  & im ){
	printf("cvPrintit not supported for double,3\n");
	return 0;
}

CVARR_TO_IMAGE_FUNC1( cvPrintit );

int main(int argc, char ** argv){
	IplImage * im = cvCreateImage( cvSize(32,32), IPL_DEPTH_64F, 3);
	IplImage * im2 = cvCreateImage( cvSize(32,32), IPL_DEPTH_64F, 1);
	cvPrintit( im );
	cvPrintit( im2 );
	float min,max;
	cvZero( im2 );
	cvTopN( im2, 50, &min, &max );
	cvTopN( im, 50, &min, &max );
	cvCopyN( im, im2 );
	printf("%f %f\n", min, max );
}
