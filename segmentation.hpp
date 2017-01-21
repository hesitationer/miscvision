#ifdef _CH_
#pragma package <opencv>
#endif

#ifndef _EiC
#include "cv.h"
#include "highgui.h"
#include <math.h>
#endif

IplImage*  image[2] = { 0, 0 }, *image0 = 0, *image1 = 0;
CvSize size;

int  w0, h0,i;
int  threshold1, threshold2;
int  l,level = 4;
int sthreshold1, sthreshold2;
int  l_comp;
int block_size = 1000;
float  parameter;
double threshold;
double rezult, min_rezult;
CvFilter filter = CV_GAUSSIAN_5x5;
CvConnectedComp *cur_comp, min_comp;
CvSeq *comp;
CvMemStorage *storage;

CvPoint pt1, pt2;
class Segmentation : public cv::Seq< CvConnectedComp > {
	CvMemStorage * m_storage;
	Segmentation() {
		m_storage = cvCreateMemStorage( 0 );
	}

	void segment(int level, int thresh1, int thresh2){
    	cvPyrSegmentation(image0, image1, m_storage, this, 
     	                 level, threshold1+1, threshold2+1);
	}

}

int main( int argc, char** argv )
{
    char* filename = argc == 2 ? argv[1] : (char*)"fruits.jpg";
    
    if( (image[0] = cvLoadImage( filename, 1)) == 0 )
        return -1;

    cvNamedWindow("Source", 0);
    cvShowImage("Source", image[0]);

    cvNamedWindow("Segmentation", 0);

    storage = cvCreateMemStorage ( block_size );

    image[0]->width &= -(1<<level);
    image[0]->height &= -(1<<level);

    image0 = cvCloneImage( image[0] );
    image1 = cvCloneImage( image[0] );
    // segmentation of the color image
    l = 1;
    threshold1 =255;
    threshold2 =30;

    ON_SEGMENT(1);

    sthreshold1 = cvCreateTrackbar("Threshold1", "Segmentation", &threshold1, 255, ON_SEGMENT);
    sthreshold2 = cvCreateTrackbar("Threshold2", "Segmentation",  &threshold2, 255, ON_SEGMENT);

    cvShowImage("Segmentation", image1);
    cvWaitKey(0);

    cvDestroyWindow("Segmentation");
    cvDestroyWindow("Source");

    cvReleaseMemStorage(&storage );

    cvReleaseImage(&image[0]);
    cvReleaseImage(&image0);
    cvReleaseImage(&image1);

    return 0;
}

#ifdef _EiC
main(1,"pyramid_segmentation.c");
#endif
