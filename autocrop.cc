#include <cv/matrix.hh>
#include <cv/cvext.h>
#include <cv/image.h>

int optimizeVert( Matrix<float> & mean, Matrix<float> & data ){
	CvMat tmplt;
	cvGetRows( &mean, &tmplt, 10, data.rows-20 );
	Matrix<float> score(data.rows-tmplt.rows+1, 1);
	cvMatchTemplate( &data, &tmplt, &score, CV_TM_CCOEFF_NORMED );

	CvPoint min,max;
	cvMinMaxLoc( &score, NULL, NULL, &min, &max );
	return max.y-10;
}

void plotHoriz( IplImage * im, Matrix<float> & vec, int y, CvScalar color){
	double min,max;
	cvMinMaxLoc(&vec, &min, &max);
	max = MAX( max+1, -min+1 );
	double scale=MIN(y, im->height-y-1)/max;
	int lastidx = cvFloor( y+vec[0][0]*scale);
	for(int i=1; i<vec.rows; i++){
		//printf("%d, %f, %f\n", x, vec[i][0]*scale, scale);
		//cvSet2D( im, i, cvFloor( x+vec[i][0]*scale ), color );
		int idx = cvFloor( y+vec[i][0]*scale );
		cvLine(im, cvPoint(i-1, lastidx), cvPoint(i, idx), color);
		lastidx = idx;
	}
}
void plotHorizLine( IplImage * im, Matrix<uchar> & vec, CvScalar color){
	for(int i=1; i<vec.rows; i++){
		if(vec[i][0] != 0){
			cvLine(im, cvPoint(i, 0), cvPoint(i, im->height ), color); 
		}
	}
} 
void plotVertLine( IplImage * im, Matrix<uchar> & vec, CvScalar color){
	for(int i=1; i<vec.rows; i++){
		if(vec[i][0] != 0){
			cvLine(im, cvPoint(0, i), cvPoint(im->width, i ), color); 
		}
	}
} 
void plotVert( IplImage * im, Matrix<float> & vec, int x, CvScalar color){
	double min,max;
	cvMinMaxLoc(&vec, &min, &max);
	max = MAX( max+1, -min+1 );
	double scale=MIN(x, im->width-x-1)/max;
	int lastidx = cvFloor( x+vec[0][0]*scale);
	for(int i=1; i<vec.rows; i++){
		//printf("%d, %f, %f\n", x, vec[i][0]*scale, scale);
		//cvSet2D( im, i, cvFloor( x+vec[i][0]*scale ), color );
		int idx = cvFloor( x+vec[i][0]*scale );
		cvLine(im, cvPoint(lastidx, i-1), cvPoint(idx, i), color);
		lastidx = idx;
	}
}
int findInRangeLeft( Matrix<float> & vec, float min, float max){
	int imin = vec.rows*0.05;
	int imax = vec.rows*0.95;
	for(int i=imin; i<imax; i++){
		if(vec[i][0]>=min && vec[i][0]<=max){
			return i;
		}
	}
	return 0;
}
int findInRangeRight( Matrix<float> & vec, float min, float max){
	int imin = vec.rows*0.05;
	int imax = vec.rows*0.95;
	for(int i=imax; i>=imin; i--){
		if(vec[i][0]>=min && vec[i][0]<=max){
			return i;
		}
	}
	return vec.rows-1;
}

template <typename T>
void sumCols( const cv::Image<T> & im, Matrix<float> & sum, CvRect roi){
	sum.realloc( roi.width, 1 );
	cvZero(&sum);
	for(int j=roi.y; j<roi.y+roi.height; j++){
		for(int i=roi.x; i<roi.x+roi.width; i++){
			sum[i][0] += im[j][i];
		}
	}
}
template <typename T>
void sumRows( const cv::Image<T> & im, Matrix<float> & sum, CvRect roi){
	sum.realloc( roi.height, 1 );
	cvZero(&sum);
	for(int j=roi.y; j<roi.y+roi.height; j++){
		for(int i=roi.x; i<roi.x+roi.width; i++){
			sum[j][0] += im[j][i];
		}
	}
}
void usage(char ** argv){
	fprintf(stderr, "USAGE: %s <src_image> <dst_image>\n", argv[0]);
}

int main(int argc, char ** argv){
	cvRedirectError( cvSegReport);
	cv::Image<uchar, 1> im;
	cv::Image<float, 1> imf, Ix, Iy;
	cv::Image<uchar, 3> im3;
	Matrix<float> vhist,hhist;
	Matrix<uchar> v_ismax, h_ismax;
	if(argc < 3){
		usage(argv);
		return -1;
	}

	assert( im.open(argv[1]) );
	im = im.resize(480, 480*im.height/im.width);
	imf = im;
	Ix.realloc( imf );
	Iy.realloc( imf );
	
	vhist.realloc( im.width, 1 );
	hhist.realloc( im.height, 1 );

	cvSobel( &imf, &Ix, 1, 0, CV_SCHARR);
	cvSobel( &imf, &Iy, 0, 1, CV_SCHARR);

	sumRows( Iy, vhist, cvRect(0, 0, im.width, im.height));
	sumCols( Ix, hhist, cvRect(0, 0, im.width, im.height));

	v_ismax.realloc( vhist.rows, vhist.cols );
	h_ismax.realloc( hhist.rows, hhist.cols );

	float min,max;
	cvMinMaxTopN(&hhist, 3, &min, &max, NULL,NULL,true); 
	printf("min=%f max=%f\n", min, max);
	int left = findInRangeLeft(hhist, min, max);
	//cvLine( &im3, cvPoint(left, 0), cvPoint(left,im3.height), CV_RGB(0,255,0));

	
	//cvInRangeS(&hhist, cvScalar(min), cvScalar(max), &h_ismax);
	//plotHorizLine( &im3, h_ismax, CV_RGB(0,255,0));

	cvMinMaxTopN(&hhist, 3, &min, &max); 
	printf("min=%f max=%f\n", min, max);
	int right = findInRangeRight(hhist, min, max);
	//cvLine( &im3, cvPoint(right, 0), cvPoint(right,im3.height), CV_RGB(0,255,0));
	//cvInRangeS(&hhist, cvScalar(min), cvScalar(max), &h_ismax);
	//plotHorizLine( &im3, h_ismax, CV_RGB(0,255,0));
	
	cvMinMaxTopN(&vhist, 5, &min, &max, NULL, NULL, true);
	printf("min=%f max=%f\n", min, max);
	int top = findInRangeLeft(vhist, min, max);
	//cvLine( &im3, cvPoint(0, top), cvPoint(im3.width, top), CV_RGB(0,255,0));
	//cvInRangeS(&vhist, cvScalar(min), cvScalar(max), &v_ismax);
	//plotVertLine( &im3, v_ismax, CV_RGB(0,255,0));
	
	cvMinMaxTopN(&vhist, 5, &min, &max);
	printf("min=%f max=%f\n", min, max );
	int bottom = findInRangeRight(vhist, min, max);
	//cvLine( &im3, cvPoint(0, bottom), cvPoint(im3.width, bottom), CV_RGB(0,255,0));
	/*
	cvInRangeS(&vhist, cvScalar(min), cvScalar(max), &v_ismax);
	plotVertLine( &im3, v_ismax, CV_RGB(0,255,0));
*/
	im3.open(argv[1]);
	left = (left * im3.width)/im.width;
	right = (right * im3.width)/im.width;
	top = (top * im3.height)/im.height;
	bottom = (bottom * im3.height)/im.height;
	
	im3.setImageROI( cvRect(left,top, right-left, bottom-top) );
	
	cv::Image<uchar, 3> dst3( right-left, bottom-top );
	cvCopy( &im3, &dst3 );
	
	dst3.save(argv[2]);
}
