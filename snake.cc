#include <stdio.h> 
 #include <math.h> 
 
 #include "cv.h" 
 #include "highgui.h" 
 #include "cvtypes.h"   // CV_TERMCRIT_ITER, M_ITER ?????????? 
 
 #define      SEGMENT            200      // ???? 
 #define      W_WIDTH            5      // ?????(WIDTH)   1, 5, 13, 25, ... 
 #define      W_HEIGHT         5      // ?????(HEIGHT)   1, 5, 13, 25, ... 
 
 IplImage *image01 = 0, *image02 = 0,  *image03 = 0;   // ????????? 
 char wndname01[] = "SNAKE??";                  // ???????? 01 
 
 int main( int argc, char** argv ) 
 { 
    CvPoint pt[SEGMENT];               // snake ??????????? 
    CvTermCriteria crit;               // ??????? ? "cvtypes.h" ?? 
    char file_name[] = "sample.jpg";      // ?????????? 
    int iteration = 10;                  // SNAKE ????? 
    CvSize win; 
    float alpha, beta, gamma; 
    int i, j; 
 
    // HighGUI ???????????? 
    cvNamedWindow( wndname01, 0 ); 
 
    // HighGUI ???????(image01)???? 
     image01 = cvLoadImage( argc > 1 ? argv[1] : file_name, 1 );   // ??? 
    image03 = cvLoadImage( argc > 1 ? argv[1] : file_name, 1 );   // ??? 
 
    // ???????????? 
    if( argc == 3 ) iteration = atoi(argv[2]); 
 
    // ??????? 
     if( !image01 ) 
     { 
         printf("Image \" %s \" was not loaded.\n", argv[1]); 
         return -1; 
     } 
 
    // HighGUI ?????????????(image02)??? 
     image02 = cvCreateImage( cvSize(image01->width, image01->height), IPL_DEPTH_8U, 1); 
     
    // ??? (image01 : BGR?? ? image02 : GRAY??) 
     cvCvtColor(image01,image02, CV_BGR2GRAY); 
     
    // Snake???????? 
    for(i=0; i<SEGMENT; i++){ 
       pt[i].x = cvRound( 0.8 * image01->width * cos(i * 6.28 / SEGMENT) / 2.0 
                                              + image01->width / 2 ); 
       pt[i].y = cvRound( 0.8 * image01->height * sin(i * 6.28 / SEGMENT) / 2.0 
                                              + image01->height / 2 ); 
    } 
 
    // ???????? 
    alpha = (float)1.3; beta = (float)0.7; gamma = (float)0.5;   // ??????? 
    //alpha = (float)1.3; beta = (float)0.1; gamma = (float)0.5; 
 
    win.width = W_WIDTH; 
    win.height = W_HEIGHT; 
    crit.type = CV_TERMCRIT_ITER;   // ??? 
    crit.max_iter = 10;            // ????? 
    //crit.maxIter = M_ITER;      // ????? 
 
    for(i=0; i<iteration; i++){ 
       // Snake??? ("cv.h" ??) 
       cvSnakeImage( image02,   // ??? (1-Channel) 
                  pt,      // ??? 
                  SEGMENT,   // ???? 
                  &alpha,   // ??? (?????????) 
                  &beta,   // ?? (??snake ??????????????) 
                  &gamma,   // ??????? 
                  1,      // int coeffUsage(????) ? DEFAULT(1) 
                  win,      // ??????????????????? 
                  crit,      // ???? 
                  1);      // ??????????????????????? DEFAULT(1) 
 
       // ????? (image01 ? image03) 
       cvCopy( image01, image03, NULL); 
 
       // ??????????? 
       for(j=0; j<SEGMENT; j++) 
       { 
          if(j < SEGMENT-1) cvLine( image03, pt[j], pt[j+1], CV_RGB(255,0,0), 2, 8 ); 
          else cvLine( image03, pt[j], pt[0], CV_RGB(255,0,0),  2, 8 ); 
       } 
 
       // ???????? 
       cvResizeWindow( wndname01, image01->width,image01->height ); 
       cvShowImage( wndname01, image03 ); 
	   cvWaitKey(-1);
    } 
 
    printf("snake ???????\n???????????\n"); 
 
    // ??????????? 
    cvWaitKey(0); 
    cvDestroyWindow(wndname01);      // ??????? 
    cvReleaseImage(&image01); 
    cvReleaseImage(&image02); 
    cvReleaseImage(&image03); 
 
    return 0; 
 } 
