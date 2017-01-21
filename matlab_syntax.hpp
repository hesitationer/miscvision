#ifndef __MATLAB_SYNTAX_H
#define __MATLAB_SYNTAX_H

#include <cv.h>
#include <highgui.h>
#include <cv/matrix.hh>

namespace matlab_syntax {


void gradient(CvArr * im, CvArr * Ix, CvArr * Iy){
	/*cvSobel(im, Ix, 1,0, 1);
	cvSobel(im, Iy, 0,1, 1);
	cvScale(Ix, Ix, 1./2);
	cvScale(Iy, Iy, 1./2);
	*/
	cvSobel(im, Ix, 1, 0, CV_SCHARR);
	cvSobel(im, Iy, 0, 1, CV_SCHARR);
	cvScale(Ix, Ix, 1./32);
	cvScale(Iy, Iy, 1./32);
	//double min,max;
	//cvMinMaxLoc(Ix, &min, &max);
	//std::cout<<"min/max "<<min<<","<<max<<std::endl;
}

template <typename data_t>
void meshgrid(float x1, float x2, float y1, float y2, Matrix<data_t> & grid){
	meshgrid(x1,1.f,x2,y1,1.f,y2,grid);
}

template<typename data_t>
void meshgrid(float x1, float xstep, float x2,
              float y1, float ystep, float y2,
			  Matrix<data_t> & result){
	int nx = cvCeil((x2-x1)/xstep);
	int ny = cvCeil((y2-y1)/ystep);
	fprintf(stdout, "%% meshgrid( %f:%f:%f, %f,%f,%f\n",x1,xstep,x2, y1,ystep,y2);
	assert(result.rows == nx*ny || fprintf(stderr, "%d,%d,%d -- %f %f %f -- %f %f %f\n", result.rows, nx, ny, x1, xstep, x2, y1, ystep, y2)==0);
	if(result.rows != nx*ny || result.cols < 2){
		if(result.cols < 2){
			result.resize(nx*ny, 2);
		}
		else{
			result.resize(nx*ny, result.cols);
		}
	}
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx; i++){
			result(j*nx+i, 0) = i*xstep+x1;
			result(j*nx+i, 1) = j*ystep+y1;
		}
	}
}

CvPoint g_ClickedPoint;
CvPoint g_MovedPoint;
bool g_WasClicked=false;
int  g_CurrentFigure = -1;
char g_CurrentName[256];
IplImage *g_CurrentImage = NULL;

char * sprintf_figure_name(char * buf, int figure){
	sprintf(buf, "Figure No. %d", figure);
	return buf;
}


void ginputCB( int event, int x, int y, int flags, void * param ){
	switch(event){
		case CV_EVENT_LBUTTONDOWN:
			if(!g_WasClicked){
				g_ClickedPoint.x = x;
				g_ClickedPoint.y = y;
				g_WasClicked=true;
			}
			break;
		case CV_EVENT_RBUTTONDOWN:
			if(!g_WasClicked){
				g_ClickedPoint.x = -1;
				g_ClickedPoint.y = -1;
				g_WasClicked=true;
			}
			break;
		case CV_EVENT_MOUSEMOVE:
			g_MovedPoint.x = x;
			g_MovedPoint.y = y;
			break;
	}
}

int figure(int figure=-1){
	char name[256];
	
	// if unspecified, find first unused
	if(figure == -1){
		for(figure=1; figure<256; figure++){
			// open new window with unused handle
			if(cvGetWindowHandle(sprintf_figure_name(name, figure))==NULL){
				cvNamedWindow(name, 0);
				break;
			}
		}
	}

	// open new window if unopen
	else if(cvGetWindowHandle(sprintf_figure_name(name, figure))==NULL){ 
		cvNamedWindow(name, 0);
	}
	
	// save some work 
	if(g_CurrentFigure == figure) return g_CurrentFigure;
	
	// unset previous mouse callbacks
	if(g_CurrentFigure > 0){
		cvSetMouseCallback(g_CurrentName, NULL);
	}

	// set active figure, name
	g_CurrentFigure = figure;
	strcpy(g_CurrentName, name);
	
	return g_CurrentFigure;
}

int gcf(){
	if(g_CurrentFigure==-1){
		return figure();
	}
	return g_CurrentFigure;
}

CvPoint ginput(){
	// install mouse handler on current figure
	cvSetMouseCallback(g_CurrentName, ginputCB);

	// reset 'clicked' indicator
	g_WasClicked=false;
	
	// wait until g_ClickNum increments
	while(! g_WasClicked ) cvWaitKey(1);

	// uninstall
	cvSetMouseCallback(g_CurrentName, NULL);
	return g_ClickedPoint;
}

void drawnow(){
	cvWaitKey(10);
}

CvPoint gmouse(){
	// install mouse handler on current figure
	cvSetMouseCallback(g_CurrentName, ginputCB);

	cvWaitKey(1);
	return g_MovedPoint;
}

// imshow
void imshow(const CvArr * im){
	gcf();
	cvShowImage(g_CurrentName, im);	
}

// imread

///
void imagesc(const CvArr * im, double min, double max){
	CvSize size = cvGetSize(im);
	
	// if im is 3-channel color -- don't scale, just do imshow
	if(cvGetElemType(im)==CV_8UC3){
		return imshow(im);
	}

	// (re) allocate temporary buffer if neccessary
	if(! g_CurrentImage ){
		g_CurrentImage = cvCreateImage(size, 8, 1);
	}
	else if(g_CurrentImage->width != size.width ||
	        g_CurrentImage->height != size.height){
		cvReleaseImage( &g_CurrentImage );
		g_CurrentImage = cvCreateImage(size, 8, 1);
	}

	// normalize image to provided values
	cvScale(im, g_CurrentImage, 255/(max-min), -min);

	// show it!
	return imshow(g_CurrentImage);
}

///  
void imagesc(const CvArr * im){
	
	double min,max;
	
	// if im is 3-channel color -- don't scale, just do imshow
	if(cvGetElemType(im)==CV_8UC3){
		return imshow(im);
	}
	
	// calculate min,max of image
	cvMinMaxLoc(im, &min, &max);
	
	return imagesc(im, min, max);
}

void pause( float delay=0 ) {
	if(delay<=0){
		cvWaitKey();
	}
	else{
		cvWaitKey( cvRound( delay*1000 ) );
	}
}

// colormap

} //namespace matlab_syntax

#endif


