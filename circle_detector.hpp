#ifndef __CV_CIRCLE_DETECTOR_HPP
#define __CV_CIRCLE_DETECTOR_HPP

#include <cv/image.h>
#include <cv/matrix.hh>
#include <cv/matrixnd.hpp>
#include <cv/contour_detector.hpp>
#include <cv/seq.hpp>
#include <cv/cvext.h>
#include <cv/cvransac.h>
#include <vector>
#include <cv/sort_by.hpp>
#include <cxmisc.h>

struct CircleVal {
	CvCircle32f * c;
	float score;
};

inline CvCircle32f cvCircleFrom3D32f( CvPoint3D32f p ){
	return cvCircle32f(p.x,p.y,p.z);
}

// 
#if 0
int calcBin(float val, float min, float max, int nbins){
}

void circleHist(cv::Image<float,1> &Ix, cv::Image<float,1> & Iy, cv::MatrixND<float, 3> & hist, 
		float xmin, float xmax, float ymin, float ymax,
		float cxmin, float cxmax, float cymin, float cymax,
		float rmin, float rmax){
	for(int y=ymin; y<ymax; y++){
		for(int x=xmin; x<xmax; x++){
			// each circle w/ center along gradient line
		}
	}
}
#endif
float sumCircle(const cv::Image<float, 1> &im, CvCircle32f circle){
	float r2 = circle.r * circle.r;
	float sum=0;
	float sz=0;
	int x0,x1,y0,y1;
	x0 = MAX( cvFloor(circle.c.x - circle.r - 1), 0);
	y0 = MAX( cvFloor(circle.c.y - circle.r - 1), 0);
	x1 = MIN( cvCeil1(circle.c.x - circle.r + 1), im.width );
	y1 = MIN( cvCeil1(circle.c.y - circle.r + 1), im.height );
	for(int j=y0; j<y1; j++){
		int xmin = cvFloor(sqrt( r2 - (j-circle.c.y)*(j-circle.c.y)));
		int xmax;
		xmax = MAX(MIN(cvCeil1(circle.c.x + xmin), im.width), 1);
		xmin = MIN(MAX(cvCeil1(circle.c.x - xmin), 0), im.width-1);
		for(int i=xmin; i<xmax; i++){
			sz++;
			sum+=im[j][i];
		}
	}
	if(sz==0) return 0;

	return sum/sz;
}

float circleError1(CvCircle32f circle, const cv::Image<float,1> &im, bool invert=false){
    float r2 = circle.r * circle.r;
    float outsum=0, insum=0;
    float outsz=0, insz=0;
	int x0,x1,y0,y1;
	x0 = MAX( cvFloor(circle.c.x - circle.r - 1), 0);
	y0 = MAX( cvFloor(circle.c.y - circle.r - 1), 0);
	x1 = MIN( cvCeil1(circle.c.x - circle.r + 1), im.width );
	y1 = MIN( cvCeil1(circle.c.y - circle.r + 1), im.height );
	for(int j=y0; j<y1; j++){
		int i;
		int xmin = cvFloor(sqrt( r2 - (j-circle.c.y)*(j-circle.c.y)));
		int xmax;
		xmax = MAX(MIN(cvCeil1(circle.c.x + xmin), im.width), 1);
		xmin = MIN(MAX(cvCeil1(circle.c.x - xmin), 0), im.width-1);
		// outside
		for(i=x0; i<xmin; i++){
			outsum += im[j][i];
			outsz ++;
		}
		// inside
		for(i=xmin; i<xmax; i++){
			insum += im[j][i];
			insz ++;
		}
		// outside again
		for(i=xmax; i<x1; i++){
			outsum += im[j][i];
			outsz ++;
		}
    }
	if(outsz<1) outsz = 1;
	if(insz<1) insz = 1;
    return (invert ? fabs( 255 - insum/insz ) + fabs( outsum/outsz ) :
			         fabs( insum/insz ) + fabs( 255 - outsum/outsz ) );
}

// compute error for circle on image im 
// assume circle is black on white unless invert==true
// error is avg difference of pixels between circle of radius c.r+1 and c.r
float circleError(CvCircle32f circle, const cv::Image<float,1> &im, bool invert=false){
	float angle_inc = 1.0 / (circle.r*M_PI*2); 
	double err=0;
	float xVec, yVec;
	int col1,row1,col2,row2;
	float v1, v2;
	int n=0;
	v1 = (invert ? 0 : 255);
	v2 = (invert ? 255 : 0);

	for (float angle = 0; angle<2*M_PI; angle+=angle_inc){
		xVec=cos(angle);
		yVec=sin(angle);
		col1 = cvFloor(circle.c.x+circle.r*xVec);
		row1 = cvFloor(circle.c.y+circle.r*yVec);
		col2 = cvFloor(circle.c.x+(circle.r+1)*xVec);
		row2 = cvFloor(circle.c.y+(circle.r+1)*yVec);
		if ((row1 >= im.height) || (row1 < 0) || (col1 < 0) || (col1 >= im.width) ||
		    (row2 >= im.height) || (row2 < 0) || (col2 < 0) || (col2 >= im.width)) continue;

		err += fabs(v1-im[row2][col2]) + fabs(v2-im[row1][col1]);
		n++;
	}

	if(n==0) n = 1024;
	return err / n ;
}

// error is avg difference of pixels between circle of radius c.r+1 and c.r
float evalCircle(CvCircle32f circle, const cv::Image<float,1> &im, bool invert=false){
	float angle_inc = 1.0 / (circle.r*M_PI*2); 
	float score=0;
	float xVec, yVec;
	int col1,row1,col2,row2;

	for (float angle = 0; angle<2*M_PI; angle+=angle_inc){
		xVec=cos(angle);
		yVec=sin(angle);
		col1 = cvFloor(circle.c.x+circle.r*xVec);
		row1 = cvFloor(circle.c.y+circle.r*yVec);
		col2 = cvFloor(circle.c.x+(circle.r+1)*xVec);
		row2 = cvFloor(circle.c.y+(circle.r+1)*yVec);
		if ((row1 >= im.height) || (row1 < 0) || (col1 < 0) || (col1 >= im.width) ||
		    (row2 >= im.height) || (row2 < 0) || (col2 < 0) || (col2 >= im.width)) continue;

		score += (im[row2][col2] - im[row1][col1]);
	}

	return score * angle_inc * (invert ? -1 : 1);
}

float evalCircle1(CvCircle32f circle, const cv::Image<float,1> &im, bool invert=false){
	// score is % of pixels for which circle of radius c.r+1
	// are greater than that of c.r
	float angle_inc = 1.0 / (circle.r*M_PI*2); 
	float score=0;
	float xVec, yVec;
	int col1,row1,col2,row2;

	for (float angle = 0; angle<2*M_PI; angle+=angle_inc){
		xVec=cos(angle);
		yVec=sin(angle);
		col1 = cvFloor(circle.c.x+circle.r*xVec);
		row1 = cvFloor(circle.c.y+circle.r*yVec);
		col2 = cvFloor(circle.c.x+(circle.r+1)*xVec);
		row2 = cvFloor(circle.c.y+(circle.r+1)*yVec);
		if ((row1 >= im.height) || (row1 < 0) || (col1 < 0) || (col1 >= im.width) ||
		    (row2 >= im.height) || (row2 < 0) || (col2 < 0) || (col2 >= im.width)) continue;

		// favors black circles on white bg by default
		score += ( (invert ? im[row1][col1] > im[row2][col2] : im[row2][col2] > im[row1][col1] ) ? 1 : 0);
	}

	return score * angle_inc;
}

// error is sum of component of gradient perpindicular to circle
float circleGradError( CvCircle32f circle, const cv::Image<float,1> &Ix, const cv::Image<float,1> &Iy) {
    //IplImage* Ix, IplImage * Iy){
    //float angle_inc = 1.0 / (circle.radius); 
    float gangle, angle_inc = 1.0 / (circle.r*M_PI*2);
    float score=0;
    float xVec, yVec;
    int col,row;
    //CvPoint2D32f gradient;
    float mag;
	int n=0;

    // 'circleness' is defined as the sum of the dot product of the
    // gradient with the directional vector along the circumference of the circle
    for (float angle = 0; angle<2*M_PI; angle+=angle_inc){
        xVec=cos(angle);
        yVec=sin(angle);
        col = cvFloor(circle.c.x+circle.r*xVec);
        row = cvFloor(circle.c.y+circle.r*yVec);
        if ((row >= Ix.height) ||
            (row < 0) || (col < 0) ||
            (col >= Ix.width)) continue;

        CvPoint2D32f gradient = cvPoint2D32f( Ix[row][col], Iy[row][col] );

        mag = sqrt(gradient.x*gradient.x + gradient.y*gradient.y);
        if(mag==0) continue;
		
		gangle = atan2( gradient.y, gradient.x );

        // calculate diff in angle
		gangle = fabs( angle - gangle );

		while(gangle > 2*M_PI) gangle-=2*M_PI;

        // sum up the dots
        // score += gangle;
		score += (gradient.x*xVec+gradient.y*yVec)/mag;
		n++;
    }
	if(n==0) return 1.0;

    // change this to negative to find white circles on black background
    return score/n; // scale by # samples 
}

float evalCircle(CvCircle32f circle, const cv::Image<float,1> &Ix, const cv::Image<float,1> &Iy){
	//IplImage* Ix, IplImage * Iy){
	//float angle_inc = 1.0 / (circle.radius); 
	float angle_inc = 1.0 / (circle.r*M_PI*2); 
	float score=0;
	float xVec, yVec;
	int col,row;
	//CvPoint2D32f gradient;
	float mag,dot;

	// 'circleness' is defined as the sum of the dot product of the
	// gradient with the directional vector along the circumference of the circle
	for (float angle = 0; angle<2*M_PI; angle+=angle_inc){
		xVec=cos(angle);
		yVec=sin(angle);
		col = cvFloor(circle.c.x+circle.r*xVec);
		row = cvFloor(circle.c.y+circle.r*yVec);
		if ((row >= Ix.height) ||
			(row < 0) || (col < 0) ||
			(col >= Ix.width)) continue;

		CvPoint2D32f gradient = cvPoint2D32f( Ix[row][col], Iy[row][col] );
		mag = sqrt(gradient.x*gradient.x + gradient.y*gradient.y);
		if(mag==0) mag=1;
		gradient.x /= mag;
		gradient.y /= mag;

		// want gradient to point in same direction as (xVec, yVec)
		// calculate the dot product
		dot = xVec*gradient.x + yVec*gradient.y;

		// sum up the dots
		score += dot;
	}

	// change this to negative to find white circles on black background
	return score*angle_inc; // scale by # samples 
}

static CvStatus icvFitCircle_32f( CvSeq * points, CvCircle32f * circle){
	CvStatus status = CV_OK;

	float c,offx,offy,sum;
	int n = points->total;
	CvSeqReader reader;
	int is_float = CV_SEQ_ELTYPE(points) == CV_32FC2;
	int i;


	offx = offy = sum = 0;
	cvStartReadSeq( points, &reader );

	// calculate center of mass
	for( i = 0; i < n; i++ )
	{
		if( !is_float )
		{
			CvPoint * p = (CvPoint *)reader.ptr;
			offx += p->x; 
			offy += p->y;
		}
		else
		{
			CvPoint2D32f * p = (CvPoint2D32f *)reader.ptr;
			offx += p->x;
			offy += p->y;
		}
		CV_NEXT_SEQ_ELEM( points->elem_size, reader );
	}

	c = 1.f / n;
	offx *= c;
	offy *= c;
	
	// calculate sum x^2 + y^2 centered at 0
	for( i = 0; i < n; i++ )
	{ 
		if( !is_float ){
			CvPoint * p = (CvPoint *)reader.ptr;
			sum += (float)((p->x-offx)*(p->x-offx) + (p->y-offy)*(p->y-offy));
		}
		else{
			CvPoint2D32f * p = (CvPoint2D32f *)reader.ptr;
			sum += (float)((p->x-offx)*(p->x-offx) + (p->y-offy)*(p->y-offy));
		}
		CV_NEXT_SEQ_ELEM( points->elem_size, reader );
	}
	sum /= n;
	*circle = cvCircle32f(offx, offy, sqrt(sum)); 

	return status;
}

int icvFitCircleRobustFitFunc( CvMat * _data, CvMat * idx, CvMat * model){
	// fit circle to first three data points
	// based on http://astronomy.swin.edu.au/~pbourke/geometry/circlefrom3/
	int * index = idx->data.i;
	Matrix<float> * data = (Matrix<float> *) _data;
	CvPoint2D32f p[3];
	int j=0;

	for(int i=0; j<3 && i<idx->rows; i++){
		if(index[i]<0) continue;
		p[j] = cvPoint2D32f( (*data)[index[i]][0], (*data)[index[i]][1] );
		assert(!isnan(p[j].x));
		assert(!isnan(p[j].y));
		assert(!isinf(p[j].x));
		assert(!isinf(p[j].y));
		j++;
	}
	assert(j==3);
	if(j!=3) return -1;

	// points are collinear, degenerate
	if(p[2].x == p[1].x && p[1].x == p[0].x){
		assert(p[2].x != p[1].x || p[1].x != p[0].x);
		return -1;
	}

	// swap order if line is vertical
	if(p[2].x == p[1].x){
		CvPoint2D32f tmp;
		CV_SWAP(p[2], p[1], tmp);
	}
	if(p[0].x == p[1].x){
		CvPoint2D32f tmp;
		CV_SWAP(p[1], p[0], tmp);
	}

	float ma = (p[2].y - p[1].y)/(p[2].x - p[1].x);
	float mb = (p[0].y - p[1].y)/(p[0].x - p[1].x);

	assert(ma!=mb);
	if(ma==mb) return -1;


	//printf("%f %f\n", ma, mb);
	//printf("%f %f\n", p[0].x, p[0].y);
	//printf("%f %f\n", p[1].x, p[1].y);
	//printf("%f %f\n", p[2].x, p[2].y);

	float cx = (ma*mb*(p[0].y-p[2].y) + mb*(p[0].x+p[1].x) - ma*(p[1].x+p[2].x))/(2*(mb-ma));
	float cy = - 1.0/ma * ( cx - (p[0].x+p[1].x)/2) + (p[0].y + p[1].y)/2;
	float r = L2(cvPoint2D32f(cx,cy), p[0]); 
	*(CvCircle32f *) model->data.fl = cvCircle32f( cx, cy, r );

	if(isnan(cx) || isnan(cy) || isnan(r) ||
	   isinf(cx) || isinf(cy) || isinf(r)){
		assert(0);
		return -1;
	}

	return 0;
}

// fit circle to given points
// data is Nx2 matrix of 2d data points
// idx is Nx1 matrix of indices into data 
// model is 3x1 matrix -- c.x,c.y,c.r
int icvFitCircleRobustFitFunc2( CvMat * data, CvMat * idx, CvMat * model){
	CvCircle32f c = cvCircle32f(0,0,0);
	int n=0;

	// calculate center of mass
	for(int i=0; i<idx->rows; i++){
		int * index = (int *)(idx->data.ptr+idx->step*i);
		if(*index<0) continue;
		float * row = (float *) cvPtr2D(data, *index, 0);
		c.c.x += row[0];
		c.c.y += row[1];
		n++;
	}
	c.c.x/=n;
	c.c.y/=n;

	// calculate radius
	for(int i=0; i<idx->rows; i++){
		int * index = (int *)(idx->data.ptr+idx->step*i);
		if(*index<0) continue;
		float * row = (float *) cvPtr2D(data, *index, 0);
		c.r += (row[0]-c.c.x)*(row[0]-c.c.x)+(row[1]-c.c.y)*(row[1]-c.c.y);
	}
	c.r /= n;
	c.r = sqrt(c.r);
	*(CvCircle32f *) model->data.fl = c;
	
	return 0;
}

// calculate quality of fit
// 
int icvFitCircleRobustDistFunc( CvMat * _data, CvMat * model, double inlierdist, CvMat * _inliers ){
    CvCircle32f c = *(CvCircle32f *) model->data.fl;
	Matrix<float> * data = (Matrix<float> *) _data;
	Matrix<float> * inliers = (Matrix<float> *) _inliers;
	int ninliers=0;

	// calculate distance from point to circle 
	for(int i=0; i<data->rows; i++){
		float dist = fabs( L2(c.c, cvPoint2D32f( (*data)[i][0], (*data)[i][1])) - c.r );
		if(dist<inlierdist){
			if(inliers) (*inliers)[i][0] = 1;
			ninliers++;
		}
		else if(inliers){
			(*inliers)[i][0]=0;
		}
	}
	return ninliers;
}

CvCircle32f cvFitCircleRobust( CvArr * array, int niterations=25, double mindist=2.0){
	CvSize size = cvGetSize( array );
	Matrix<int> inliers(size.height, 1);
	Matrix<float> model(3,1);

	//cvRansac( array, &inliers, &model, icvFitCircleRobustFitFunc, icvFitCircleRobustDistFunc, NULL, niterations, 3, mindist, niterations);
	return *(CvCircle32f *) model[0];
}

CvCircle32f cvFitCircle(CvArr * array){
	CvCircle32f circle; // x,y,radius

	CV_FUNCNAME( "cvFitCircle" );

	memset( &circle, 0, sizeof(circle));

	__BEGIN__;

	CvContour contour_header;
	CvSeq* ptseq = 0;
	CvSeqBlock block;

	if( CV_IS_SEQ( array ))
	{
		ptseq = (CvSeq*)array;
		if( !CV_IS_SEQ_POINT_SET( ptseq ))
			CV_ERROR( CV_StsBadArg, "Unsupported sequence type" );
	}
	else
	{
		CV_CALL( ptseq = cvPointSeqFromMat(
					CV_SEQ_KIND_GENERIC, array, &contour_header, &block ));
	}

	if( ptseq->total < 6 )
		CV_ERROR( CV_StsBadSize, "Number of points should be >= 6" );

	CV_CALL( icvFitCircle_32f( ptseq, &circle ));

	__END__;

	return circle;
}

int cmp_circle_val(const void* _a, const void* _b, void* _userdata){
	CircleVal*a=(CircleVal*)_a,*b=(CircleVal*)_b;
	return cvFloor(a->score - b->score);  // lower score better
}

int cmp_circle_area(const void* _a, const void* _b, void* _userdata){
	CvCircle32f *a=(CvCircle32f*)_a,*b=(CvCircle32f*)_b;
	return cvFloor(-(a->r*a->r - b->r*b->r));
}

float cvCirclePointDist(CvCircle32f c, CvSeq * seq){
	CvSeqReader reader;
	cvStartReadSeq( seq, &reader );
	int n = seq->total;
	float sum = 0;
	for(int i=0; i<n; i++){
		CvPoint * p = (CvPoint*) reader.ptr;
		sum += fabs( L2(c.c, cvPointTo32f(*p)) - c.r );
		CV_NEXT_SEQ_ELEM( seq->elem_size, reader );
	}
	return sum/n;
}
CvPoint2D32f avgPoint( CvSeq * seq ){
	CvSeqReader reader;
	cvStartReadSeq( seq, &reader );
	CvPoint2D32f avg = cvPoint2D32f(0,0);
	for(int i=0; i<seq->total; i++){
		avg.x += ((CvPoint*) reader.ptr)->x;
		avg.y += ((CvPoint*) reader.ptr)->y;
		CV_NEXT_SEQ_ELEM( seq->elem_size, reader );
	}
	avg.x /= seq->total;
	avg.y /= seq->total;
	return avg;
}

class CircleDetector : public ContourDetector {
protected:
	std::vector<CvCircle32f> m_circles;
	//CvMemStorage * m_storage;
public:
	typedef std::vector<CvCircle32f>::iterator iterator;
	typedef std::vector<CvCircle32f>::const_iterator const_iterator;

	CircleDetector():ContourDetector(){
		//m_storage = cvCreateMemStorage(0);
		//m_circles.init(m_storage);
	}
	void detect(IplImage * im, double thresh, bool binary=false){
		std::vector<float> scores;
		CvCircle32f c;
		//CvPoint* PointArray;
		//CvPoint2D32f* PointArray2D32f;
		// first clear old circles
		this->clear();

		ContourDetector::detect(im, thresh, binary);
		m_circles.clear();

		// This approximate each contour by an circle.
		for(CvSeq* cont=m_contours;cont;cont = cont->h_next)
		{   
			int count = cont->total; // This is number point in contour
			float score;
			// Number of points must be more than or equal to 6 (for cvFitCircle_32f).        
			if( count < 6 ){
				c.c = avgPoint( cont );
				c.r = (double) count / 12;  // want maximum radius for this case to be .5
				score = L2( cvPoint2D32f(0,0), cvPoint2D32f(im->width, im->height) ); 
			}
			else{
				// Fits circle to current contour.
				c = cvFitCircle(cont); //PointArray2D32f);
				// score = -count; // M_PI*cval.c.r*cval.c.r; //cvCirclePointDist(cval.c, cont);
				// if(c.r<1) c.r = 1;
				score = cvCirclePointDist(c, cont);
			}

			//printf("c.c=%f,%f c.r=%f score=%f\n", c.c.x, c.c.y, c.r, score);
			m_circles.push_back( c );
			scores.push_back( score );
		}
		sort_by( m_circles.begin(), m_circles.end(), scores.begin(), scores.end() );
	}
	const_iterator begin() const{
		return m_circles.begin();
	}
	const_iterator end() const {
		return m_circles.end();
	}
	iterator begin(){
		return m_circles.begin();
	}
	iterator end(){
		return m_circles.end();
	}
	static CvCircle32f refine(CvCircle32f circle, IplImage * Ix, IplImage * Iy, CvSize win=cvSize(5,5), float rwin=5){
		assert(Ix->depth==IPL_DEPTH_32F && Iy->depth==IPL_DEPTH_32F && 
			   Ix->nChannels==1 && Iy->nChannels==1);
		return CircleDetector::refine(circle, *(cv::Image<float,1> *)Ix, *(cv::Image<float,1> *)Iy, win, rwin);
	}

	static CvCircle32f refine(CvCircle32f circle, const cv::Image<float,1> & Ix, const cv::Image<float,1>& Iy, CvSize win=cvSize(5,5), float rwin=5, bool invert=false){
		CvCircle32f result;
		float val;
		CircleDetector::refine(circle, Ix, Iy, result, val, win, rwin, invert);
		return result;
	}

	static CvCircle32f refine(CvCircle32f circle, const cv::Image<float, 1> & im, CvSize win=cvSize(5,5), float rwin=5, bool invert=false){
		CvCircle32f result;
		float val;
		CircleDetector::refine(circle, im, result, val, win, rwin, invert);
		return result;
	}
	static CvCircle32f subpixel_refine(CvCircle32f circle, const cv::Image<float, 1> & im, bool invert=false){
		int dims[] = {3, 3, 3};
        cv::MatrixND< float, 3 > lut( dims );
		CvCircle32f c;
		if(circle.c.x<1 || circle.c.x>=im.width-1 || circle.c.y<1 || circle.c.y>=im.height-1 || circle.r<2) return circle;
        for(int y=-1; y<=1; y++){
            for(int x=-1; x<=1; x++){
                for(int r=-1; r<=1; r++){
                    c.c.x = circle.c.x + x;
                    c.c.y = circle.c.y + y;
                    c.r = circle.r + r;
                    lut(y+1,x+1,r+1) = circleError(c, im, invert);
                }
            }
        }
		c = subpixel_refine( lut );
		if(fabs(c.c.x)<=.5 && fabs(c.c.y)<=.5 && fabs(c.r) <=.5){
			return cvCircle32f( circle.c.x + c.c.x, circle.c.y + c.c.y, circle.r + c.r );
		}
		else return circle;
	}

	static CvCircle32f subpixel_refine(CvCircle32f circle, const cv::Image<float, 1> & Ix, const cv::Image<float,1> & Iy, bool invert=false){
		int dims[] = {3, 3, 3};
        cv::MatrixND< float, 3 > lut( dims );
		CvCircle32f c;
		if(circle.c.x<1 || circle.c.x>=Ix.width-1 || circle.c.y<1 || circle.c.y>=Ix.height-1 || circle.r<2) return circle;
        for(int y=-1; y<=1; y++){
            for(int x=-1; x<=1; x++){
                for(int r=-1; r<=1; r++){
                    c.c.x = circle.c.x + x;
                    c.c.y = circle.c.y + y;
                    c.r = circle.r + r;
                    lut(y+1,x+1,r+1) = circleGradError(c, Ix, Iy);
					if(invert) lut(y+1,x+1,r+1) = 1 - lut(y+1,x+1,r+1);
                }
            }
        }
		c = subpixel_refine( lut );
		if(fabs(c.c.x)<=.5 && fabs(c.c.y)<=.5 && fabs(c.r) <=.5){
			return cvCircle32f( circle.c.x + c.c.x, circle.c.y + c.c.y, circle.r + c.r );
		}
		else return circle;
	}

	// refine location of circle using hough-like optimization 
	static void refine(CvCircle32f circle, const cv::Image<float,1> & im, CvCircle32f & result, 
			float & value, CvSize win, float rwin, bool invert){
		int maxR = cvFloor(circle.r + rwin);
		int minR = MAX(1, cvFloor(circle.r - rwin));
		int minX = MAX(0, cvFloor(circle.c.x-win.width));
		int maxX = cvFloor(MIN(im.width-1, circle.c.x+win.width));
		int minY = MAX(0, cvFloor(circle.c.y-win.height));
		int maxY = cvFloor(MIN(im.height-1, circle.c.y+win.height));
		float maxVal = circleError(circle, im, invert);
		CvCircle32f maxC=circle;
		CvCircle32f c;
		for(c.c.y = minY; c.c.y <= maxY; c.c.y++){
			for(c.c.x = minX; c.c.x <= maxX; c.c.x++){
				for(c.r=minR; c.r<=maxR; c.r++){
					float v = circleError(c, im, invert);
					if(!isinf(v) && !isnan(v) && v < maxVal){
						maxVal = v;
						maxC = c;
					}
				}
			}
		}
		result = maxC;
		value = maxVal;
	}

	// refine location of circle using hough-like optimization with a gradient scoring function
	static void refine(CvCircle32f circle, const cv::Image<float,1> & Ix, const cv::Image<float,1> &Iy, CvCircle32f & result, float & value, CvSize win, float rwin, bool invert){
		float scale= (invert ? -1 : 1);
		int maxR = cvFloor(circle.r + rwin);
		int minR = MAX(1, cvFloor(circle.r - rwin));
		int minX = MAX(0, cvFloor(circle.c.x-win.width));
		int maxX = cvFloor(MIN(Ix.width-1, circle.c.x+win.width));
		int minY = MAX(0, cvFloor(circle.c.y-win.height));
		int maxY = cvFloor(MIN(Ix.height-1, circle.c.y+win.height));
		float maxVal = circleGradError(circle, Ix, Iy)*scale;
		CvCircle32f maxC=circle;
		CvCircle32f c;
		for(c.c.y = minY; c.c.y <= maxY; c.c.y++){
			for(c.c.x = minX; c.c.x <= maxX; c.c.x++){
				for(c.r=minR; c.r<=maxR; c.r++){
					float v = circleGradError(c, Ix, Iy) * scale;
					if(!isinf(v) && !isnan(v) && v < maxVal){
						maxVal = v;
						maxC = c;
					}
				}
			}
		}

		result = maxC;
		value = maxVal;
	}
	
	/// data is a 3x3x3 array where 1,1,1 is the current sample point
	static CvCircle32f subpixel_refine( cv::MatrixND<float, 3> & data ) {
		float dx,dy,dr,dxx,dyy,drr,dxy,dxr,dyr;

		// first derivatives df/dx ~=~ f(x+1) - f(x-1)
		dx = data(1,2,1) - data(1,0,1);
		dy = data(2,1,1) - data(0,1,1);
		dr = data(1,1,2) - data(1,1,0);
		
		// second derivatives d^2f/dx^2 ~=~ f(x-1) - 2*f(x) + f(x+1)
		dxx = data(1,0,1) - 2*data(1,1,1) + data(1,2,1);
		dyy = data(0,1,1) - 2*data(1,1,1) + data(2,1,1);
		drr = data(1,1,0) - 2*data(1,1,1) + data(1,1,2);
		
		// second partials df/dxy ~= .25( f(x+1,y+1) - f(x-1, y+1) - f(x+1, y-1) + f(x-1,y-1) )
		dxy = .25 * (data(2,2,1) + data(0,0,1) - data(0,2,1) - data(2,0,1));
		dxr = .25 * (data(1,2,2) + data(1,0,0) - data(1,2,0) - data(1,0,2));
		dyr = .25 * (data(2,1,2) + data(0,1,0) - data(2,1,0) - data(0,1,2));

		return subpixel_refine( dx, dy, dr, dxx, dyy, drr, dxy, dxr, dyr);
	}

	static CvCircle32f subpixel_refine(float dx, float dy, float dr, float dxx, float dyy, float drr, float dxy, float dxr, float dyr){
		Matrix<float> H(3,3);
		H(0,0) = dxx; H(1,1) = dyy; H(2,2) = drr;
		H(0,1) = H(1,0) = dxy;
		H(0,2) = H(2,0) = dxr;
		H(1,2) = H(2,1) = dyr;

		Matrix<float> d(3,1);
		d(0,0) = -dx;
		d(1,0) = -dy;
		d(2,0) = -dr;

		d = H.inv()*d;
		
		printf("dp = [ %f %f %f ] \n", d(0,0), d(1,0), d(2,0) );
		return cvCircle32f( d(0,0), d(1,0), d(2,0) );
	}
	
	// refine location of a pair of circles assumed to have same radius, but different locations
	static void refine_pair(CvCircle32f c1, CvCircle32f c2, const cv::Image<float,1> & im, CvCircle32f & res1, CvCircle32f &res2, float & val, CvSize win=cvSize(5,5), float rwin=5, bool invert=false){
		int maxR = cvFloor(MAX(c1.r,c2.r) + rwin);
		int minR = MAX(1, cvFloor(MIN(c1.r,c2.r) - rwin));
		
		float maxVal = 0;//(evalCircle(c1, Ix, Iy) * evalCircle(c2, Ix, Iy));
		float v1,v2,v;
		CvCircle32f maxc1=c1;
		CvCircle32f maxc2=c2;
		CvCircle32f tc1,tc2;
		for(int r=minR; r<maxR; r++){
			c1.r = r;
			c2.r = r;
			CircleDetector::refine( c1, im, tc1, v1, win, 0, invert ); 
			CircleDetector::refine( c2, im, tc2, v2, win, 0, invert );
			v = (v1+1)*(v2+1);
			if(r==minR || (!isinf(v) && !isnan(v) && v < maxVal)){
				maxVal = v;
				maxc1 = tc1;
				maxc2 = tc2;
			}
		}
		res1 = maxc1;
		res2 = maxc2;
		val = -maxVal;
	}

	// refine location of a pair of circles assumed to have same radius, but different locations
	static void refine_pair(CvCircle32f c1, CvCircle32f c2, const cv::Image<float,1> & Ix, const cv::Image<float,1>& Iy, CvCircle32f & res1, CvCircle32f &res2, float & val, CvSize win=cvSize(5,5), float rwin=5, bool invert=false){
		int maxR = cvFloor(MAX(c1.r,c2.r) + rwin);
		int minR = MAX(1, cvFloor(MIN(c1.r,c2.r) - rwin));
		
		float maxVal = 0;//(evalCircle(c1, Ix, Iy) * evalCircle(c2, Ix, Iy));
		float v1,v2,v;
		CvCircle32f maxc1=c1;
		CvCircle32f maxc2=c2;
		CvCircle32f tc1,tc2;
		for(int r=minR; r<maxR; r++){
			c1.r = r;
			c2.r = r;
			CircleDetector::refine( c1, Ix, Iy, tc1, v1, win, 0, invert ); 
			CircleDetector::refine( c2, Ix, Iy, tc2, v2, win, 0, invert );
			v = v1*v2;
			if(!isinf(v) && !isnan(v) && v > maxVal){
				maxVal = v;
				maxc1 = tc1;
				maxc2 = tc2;
			}
		}
		res1 = maxc1;
		res2 = maxc2;
		val = -maxVal;
	}
	void draw(IplImage * im, CvScalar color){
		for(size_t i=0; i<m_circles.size(); i++){
			CvCircle32f c = m_circles[i];
			cvCircle(im, cvPointFrom32f(c.c), cvCeil1(c.r), color, 1, CV_AA);
		}
	}

	CvCircle32f getCircle(int i){
		return m_circles[i];
		//(CvPoint3D32f *) cvGetSeqElem(m_circles, i);
	}
	float getScore(int i){
		return 0;
		//return m_circles[i]core;
	}
	CvCircle32f & operator[](int i){
		return m_circles[i];
	}
	int getNumCircles(){
		return m_circles.size();//->total;
	}
	int size(){
		return m_circles.size();
	}
	void clear(){
		ContourDetector::clear();
		//cvClearSeq(m_circles);
		m_circles.clear();
	}
};

#endif //__CV_CIRCLE_DETECTOR_HPP
