#include <stdio.h>
#include <highgui.h>
#include <cv/cvext.h>

CvPoint points[2];
int npoints=0;

void cb( int event, int x, int y, int flags, void * param ){
	if(event==CV_EVENT_LBUTTONUP){
		if(npoints<2){
			points[npoints] = cvPoint(x,y);
			npoints++;
		}
	}
}

void cvRotateScale( CvArr * src, CvArr * dst, CvPoint2D32f OrigA, CvPoint2D32f OrigB, CvPoint2D32f DestA, CvPoint2D32f DestB ){
	float scale;
	float theta;
	float map[6];
	float cosTheta,sinTheta;
	CvMat src_mat,dst_mat,map_mat;
	CvSize sz = cvGetSize(dst);
	CvPoint2D32f OrigMid, DestMid;

	map_mat = cvMat( 2, 3, CV_32F, map );

	scale = sqrt((OrigB.y - OrigA.y)*(OrigB.y - OrigA.y)+(OrigB.x-OrigA.x)*(OrigB.x-OrigA.x)) /
		 sqrt((DestB.y - DestA.y)*(DestB.y - DestA.y)+(DestB.x-DestA.x)*(DestB.x-DestA.x));

	OrigMid = cvPoint2D32f( OrigA.x + .5*(OrigB.x-OrigA.x), OrigA.y + .5*(OrigB.y-OrigA.y) );
	DestMid = cvPoint2D32f( DestA.x + .5*(DestB.x-DestA.x), DestA.y + .5*(DestB.y-DestA.y) );

	theta = atan2(OrigB.y - OrigA.y, OrigB.x - OrigA.x) - atan2(DestB.y - DestA.y, DestB.x - DestA.x);
	cosTheta = cos(theta);
	sinTheta = sin(theta);

	map[0] = map[4] = cosTheta*scale;
	map[1] = -sinTheta*scale;
	map[3] = sinTheta*scale;

	map[2] = 0;
	map[5] = 0;
	map[2] = OrigMid.x - map[0]*(DestMid.x) - map[1]*(DestMid.y);
	map[5] = OrigMid.y - map[3]*(DestMid.x) - map[4]*(DestMid.y);

	// x' = x-(width(dst)-1)*0.5
	//     // y' = y-(height(dst)-1)*0.5
	//         // dst(x,y) = src( A11x' + A12y' + b1, A21x'+A22y' + b2)
//	cvGetQuadrangleSubPix( src, dst, &map_mat );
	cvWarpAffine( src, dst, &map_mat, CV_WARP_INVERSE_MAP+CV_INTER_LINEAR );
}

#define DIST(p,q) (sqrt( (p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y)))

int main(int argc, char ** argv){
	cvNamedWindow("image", 1);
	cvNamedWindow("rotated", 0);
	IplImage * im = cvLoadImage( argv[1], 1 );
	IplImage * clone  = cvCloneImage( im );
	IplImage * rotated = cvCreateImage(cvSize(200,200), 8, 3); 

	cvSetMouseCallback("image", cb);
	while(1){
		if(npoints==1){
			cvCopy(im, clone);
			cvCircle( clone, points[0], 1, CV_RGB(0,255,0) );
		}
		if(npoints==2){
			float d = DIST( points[0], points[1] );
			float angle = atan2( points[1].y-points[0].y, points[1].x-points[0].x );
			float scale = rotated->width*.5/d;
			printf("scale=%f angle=%f\n", scale, angle/M_PI*180);
			cvRotateScale( im, rotated, cvPointTo32f(points[0]), cvPointTo32f(points[1]), cvPoint2D32f(rotated->width/4, rotated->height*3/4), cvPoint2D32f(rotated->width*3/4, rotated->height*3/4) );
			cvShowImage("rotated", rotated);
			cvCircle( clone, points[1], 1, CV_RGB(0,255,0) );
			cvLine( clone, points[0], points[1], CV_RGB(0,255,0) );
			npoints=0;
		}
		cvShowImage("image", clone);
		if((char)cvWaitKey(100)=='q') break;
	}
}
