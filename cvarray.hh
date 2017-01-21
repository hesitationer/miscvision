template <typename data_t> class Array : CvArr {
	virtual ~CvArr = 0;
CV_INLINE  void  cvDecRefData( CvArr* arr );
CV_INLINE  void  cvDecRefData( CvArr* arr )
CV_INLINE  int  cvIncRefData( CvArr* arr );
CV_INLINE  int  cvIncRefData( CvArr* arr )
OPENCVAPI  int cvGetElemType( const CvArr* arr );
OPENCVAPI  int cvGetDims( const CvArr* arr, int* sizes CV_DEFAULT(NULL) );
OPENCVAPI  int cvGetDimSize( const CvArr* arr, int index );
OPENCVAPI uchar* cvPtr1D( const CvArr* arr, int idx1, int* type CV_DEFAULT(NULL));
OPENCVAPI uchar* cvPtr2D( const CvArr* arr, int idx1, int idx2, int* type CV_DEFAULT(NULL) );
OPENCVAPI uchar* cvPtr3D( const CvArr* arr, int idx1, int idx2, int idx3,
OPENCVAPI uchar* cvPtrND( const CvArr* arr, int* idx, int* type CV_DEFAULT(NULL) );
OPENCVAPI CvScalar cvGet1D( const CvArr* arr, int idx1 );
OPENCVAPI CvScalar cvGet2D( const CvArr* arr, int idx1, int idx2 );
OPENCVAPI CvScalar cvGet3D( const CvArr* arr, int idx1, int idx2, int idx3 );
OPENCVAPI CvScalar cvGetND( const CvArr* arr, int* idx );
OPENCVAPI double cvGetReal1D( const CvArr* arr, int idx1 );
OPENCVAPI double cvGetReal2D( const CvArr* arr, int idx1, int idx2 );
OPENCVAPI double cvGetReal3D( const CvArr* arr, int idx1, int idx2, int idx3 );
OPENCVAPI double cvGetRealND( const CvArr* arr, int* idx );
OPENCVAPI void cvSet1D( CvArr* arr, int idx1, CvScalar value );
OPENCVAPI void cvSet2D( CvArr* arr, int idx1, int idx2, CvScalar value );
OPENCVAPI void cvSet3D( CvArr* arr, int idx1, int idx2, int idx3, CvScalar value );
OPENCVAPI void cvSetND( CvArr* arr, int* idx, CvScalar value );
OPENCVAPI void cvSetReal1D( CvArr* arr, int idx1, double value );
OPENCVAPI void cvSetReal2D( CvArr* arr, int idx1, int idx2, double value );
OPENCVAPI void cvSetReal3D( CvArr* arr, int idx1,
OPENCVAPI void cvSetRealND( CvArr* arr, int* idx, double value );
OPENCVAPI void cvClearND( CvArr* arr, int* idx );
OPENCVAPI CvMat* cvGetMat( const CvArr* src, CvMat* header,
OPENCVAPI IplImage* cvGetImage( const CvArr* arr, IplImage* img );
OPENCVAPI CvArr* cvReshapeMatND( const CvArr* arr,
                                 int sizeof_header, CvArr* header,
OPENCVAPI CvMat* cvReshape( const CvArr* arr, CvMat* header,
OPENCVAPI void cvRepeat( const CvArr* src, CvArr* dst );
OPENCVAPI void  cvCreateData( CvArr* arr );
OPENCVAPI void  cvReleaseData( CvArr* arr );
OPENCVAPI void  cvSetData( CvArr* arr, void* data, int step );
OPENCVAPI void cvGetRawData( const CvArr* arr, uchar** data,
OPENCVAPI  CvSize cvGetSize( const CvArr* arr );
OPENCVAPI  void  cvCopy( const CvArr* src, CvArr* dst,
                         const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvSet( CvArr* arr, CvScalar scalar,
                        const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvSetZero( CvArr* mat );
OPENCVAPI  void  cvConvertScale( const CvArr *src, CvArr *dst,
OPENCVAPI  void  cvAdd( const CvArr* srcA, const CvArr* srcB, CvArr* dst,
                        const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI  void  cvAddS( const CvArr* src, CvScalar value, CvArr* dst,
                         const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI  void  cvSub( const CvArr* srcA, const CvArr* srcB, CvArr* dst,
                        const CvArr* mask CV_DEFAULT(NULL));
CV_INLINE  void  cvSubS( const CvArr* src, CvScalar value, CvArr* dst,
                         const CvArr* mask CV_DEFAULT(NULL));
CV_INLINE  void  cvSubS( const CvArr* src, CvScalar value, CvArr* dst,
                         const CvArr* mask )
OPENCVAPI  void  cvSubRS( const CvArr* src, CvScalar value, CvArr* dst,
                          const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI  void  cvMul( const CvArr* srcA, const CvArr* srcB,
                        CvArr* dst, double scale CV_DEFAULT(1) );
OPENCVAPI  void  cvDiv( const CvArr* srcA, const CvArr* srcB,
                        CvArr* dst, double scale CV_DEFAULT(1));
OPENCVAPI  void  cvScaleAdd( const CvArr* srcA, CvScalar scale,
                             const CvArr* srcB, CvArr* dst );
OPENCVAPI  void  cvAddWeighted( const CvArr* srcA, double alpha,
                                const CvArr* srcB, double beta,
                                double gamma, CvArr* dst );
OPENCVAPI  double  cvDotProduct( const CvArr* srcA, const CvArr* srcB );
OPENCVAPI void cvAnd( const CvArr* src1, const CvArr* src2,
                      CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvAndS( const CvArr* src, CvScalar value,
                       CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvOr( const CvArr* src1, const CvArr* src2,
                     CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvOrS( const CvArr* src, CvScalar value,
                      CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvXor( const CvArr* src1, const CvArr* src2,
                      CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvXorS( const CvArr* src, CvScalar value,
                       CvArr* dst, const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI void cvNot( const CvArr* src, CvArr* dst );
OPENCVAPI void cvInRange( const CvArr* src, const CvArr* lower,
                          const CvArr* upper, CvArr* dst );
OPENCVAPI void cvInRangeS( const CvArr* src, CvScalar lower,
                           CvScalar upper, CvArr* dst );
OPENCVAPI void cvCmp( const CvArr* src1, const CvArr* src2, CvArr* dst, int cmpOp );
OPENCVAPI void cvCmpS( const CvArr* src1, double scalar, CvArr* dst, int cmpOp );
OPENCVAPI void cvMin( const CvArr* src1, const CvArr* src2, CvArr* dst );
OPENCVAPI void cvMax( const CvArr* src1, const CvArr* src2, CvArr* dst );
OPENCVAPI void cvMinS( const CvArr* src, double scalar, CvArr* dst );
OPENCVAPI void cvMaxS( const CvArr* src, double scalar, CvArr* dst );
OPENCVAPI  void  cvCartToPolar( const CvArr* x, const CvArr* y,
                                CvArr* magnitude, CvArr* angle CV_DEFAULT(NULL),
OPENCVAPI  void  cvPolarToCart( const CvArr* magnitude, const CvArr* angle,
                                CvArr* x, CvArr* y,
OPENCVAPI  void  cvPow( const CvArr* src, CvArr* dst, double power );
OPENCVAPI  void  cvExp( const CvArr* src, CvArr* dst );
OPENCVAPI  void  cvLog( const CvArr* src, CvArr* dst );
OPENCVAPI  int  cvCheckArr( const CvArr* arr, int flags CV_DEFAULT(0),
OPENCVAPI  void  cvRand( CvRandState* state, CvArr* arr );
OPENCVAPI  void  cvCrossProduct( const CvArr* srcA, const CvArr* srcB, CvArr* dst );
OPENCVAPI  void  cvMatMulAdd( const CvArr* srcA, const CvArr* srcB,
                              const CvArr* srcC, CvArr* dst );
OPENCVAPI  void  cvGEMM( const CvArr* srcA, const CvArr* srcB, double alpha,
                         const CvArr* srcC, double beta, CvArr* dst,
OPENCVAPI  void  cvMatMulAddS( const CvArr* src, CvArr* dst,
OPENCVAPI void cvMulTransposed( const CvArr* srcarr,
                                CvArr* dstarr, int order );
OPENCVAPI  void  cvTranspose( const CvArr* src, CvArr* dst );
OPENCVAPI  void  cvFlip( const CvArr* src, CvArr* dst CV_DEFAULT(NULL),
OPENCVAPI  void   cvSVD( CvArr* A, CvArr* W CV_DEFAULT(NULL),
                         CvArr* U CV_DEFAULT(NULL),
                         CvArr* V CV_DEFAULT(NULL),
OPENCVAPI  void   cvSVBkSb( const CvArr* warr, const CvArr* uarr,
                            const CvArr* varr, const CvArr* barr,
                            CvArr* xarr, int flags );
OPENCVAPI  double  cvInvert( const CvArr* src, CvArr* dst,
OPENCVAPI  int  cvSolve( const CvArr* A, const CvArr* b, CvArr* x,
OPENCVAPI  double cvDet( const CvArr* mat );
OPENCVAPI  CvScalar cvTrace( const CvArr* mat );
OPENCVAPI  void  cvEigenVV( CvArr* src, CvArr* evects, CvArr* evals, double eps );
OPENCVAPI  void  cvSetIdentity( CvArr* mat, CvScalar value CV_DEFAULT(cvScalar(1)) );
OPENCVAPI  void  cvPerspectiveTransform( const CvArr* src, CvArr* dst, const CvArr* mat );
OPENCVAPI  void  cvCalcCovarMatrix( const CvArr** vects, CvArr* covarMatrix, CvArr* avg );
OPENCVAPI  double  cvMahalanobis( const CvArr* srcA, const CvArr* srcB, CvArr* mat );
OPENCVAPI  CvScalar  cvSum( const CvArr* arr );
OPENCVAPI  int  cvCountNonZero( const CvArr* arr );
OPENCVAPI  CvScalar  cvAvg( const CvArr* arr, const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvAvgSdv( const CvArr* arr, CvScalar* mean, CvScalar* std_dev,
                           const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvMinMaxLoc( const CvArr* arr, double* min_val, double* max_val,
                              const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  double  cvNorm( const CvArr* arrA, const CvArr* arrB CV_DEFAULT(NULL),
                           const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvDFT( const CvArr* src, CvArr* dst, int flags );
OPENCVAPI  void  cvMulCcs( const CvArr* srcA, const CvArr* srcB, CvArr* dst );
OPENCVAPI  void  cvDCT( const CvArr* src, CvArr* dst, int flags );
OPENCVAPI  void*  cvCvtSeqToArray( CvSeq* seq, CvArr* arr,
OPENCVAPI  void  cvSeqInsertSlice( CvSeq* seq, int index, const CvArr* from_arr );
OPENCVAPI  void cvLUT( const CvArr* srcarr, CvArr* dstarr, const CvArr* lutarr );
OPENCVAPI  void cvSmooth( const CvArr* srcarr, CvArr* dstarr,
OPENCVAPI void cvIntegral( const CvArr* image, CvArr* sumImage,
                           CvArr* sumSqImage CV_DEFAULT(NULL),
                           CvArr* tiltedSumImage CV_DEFAULT(NULL));
OPENCVAPI  void  cvPyrDown( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvPyrUp( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvCalcPyramid( const CvArr* src, CvArr* container,
OPENCVAPI void cvSobel( const CvArr* src, CvArr* dst,
OPENCVAPI void cvLaplace( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvCvtColor( const CvArr* src, CvArr* dst, int colorCvtCode );
OPENCVAPI  void  cvResize( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvErode( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvDilate( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvMorphologyEx( const CvArr* src, CvArr* dst,
                                 CvArr* temp, IplConvKernel* element,
OPENCVAPI void cvMoments( const CvArr* arr, CvMoments* moments, int binary CV_DEFAULT(0));
OPENCVAPI  void  cvLine( CvArr* arr, CvPoint pt1, CvPoint pt2,
OPENCVAPI  void  cvLineAA( CvArr* arr, CvPoint pt1, CvPoint pt2,
OPENCVAPI  void  cvRectangle( CvArr* arr, CvPoint pt1, CvPoint pt2,
OPENCVAPI  void  cvCircle( CvArr* arr, CvPoint center, int radius,
OPENCVAPI  void  cvCircleAA( CvArr* arr, CvPoint center, int radius,
OPENCVAPI  void  cvEllipse( CvArr* arr, CvPoint center, CvSize axes,
CV_INLINE  void  cvEllipseBox( CvArr* arr, CvBox2D box,
CV_INLINE  void  cvEllipseBox( CvArr* arr, CvBox2D box,
OPENCVAPI  void  cvEllipseAA( CvArr* arr, CvPoint center, CvSize axes,
OPENCVAPI  void  cvFillConvexPoly( CvArr* arr, CvPoint* pts, int npts, double color );
OPENCVAPI  void  cvFillPoly( CvArr* arr, CvPoint** pts,
OPENCVAPI  void  cvPolyLine( CvArr* arr, CvPoint** pts, int* npts, int contours,
OPENCVAPI  void  cvPolyLineAA( CvArr* arr, CvPoint** pts, int* npts, int contours,
OPENCVAPI  void  cvPutText( CvArr* arr, const char* text, CvPoint org,
OPENCVAPI  int  cvInitLineIterator( const CvArr* arr, CvPoint pt1, CvPoint pt2,
OPENCVAPI  int  cvSampleLine( const CvArr* arr, CvPoint pt1, CvPoint pt2, void* buffer,
OPENCVAPI  void  cvGetRectSubPix( const CvArr* src, CvArr* dst, CvPoint2D32f center );
OPENCVAPI  void  cvGetQuadrangleSubPix( const CvArr* src, CvArr* dstarr,
                                        const CvArr* matrixarr,
OPENCVAPI  void  cvMatchTemplate( const CvArr* arr, const CvArr* templ,
                                  CvArr* result, int method );
OPENCVAPI  float  cvCalcEMD2( const CvArr* signature1,
                              const CvArr* signature2,
                              const CvArr* cost_matrix CV_DEFAULT(0),
                              CvArr* flow CV_DEFAULT(0),
OPENCVAPI  int  cvFindContours( CvArr* arr, CvMemStorage* storage,
OPENCVAPI  CvContourScanner   cvStartFindContours( CvArr* arr, CvMemStorage* storage,
OPENCVAPI void  cvDrawContours( CvArr *img, CvSeq* contour,
OPENCVAPI  void  cvAbsDiff( const CvArr* srcA, const CvArr* srcB, CvArr* dst );
OPENCVAPI  void  cvAbsDiffS( const CvArr* src, CvArr* dst, CvScalar value );
OPENCVAPI  void  cvCalcOpticalFlowLK( const CvArr* srcA, const CvArr* srcB,
                                      CvSize winSize, CvArr* velx, CvArr* vely );
OPENCVAPI  void  cvCalcOpticalFlowBM( const CvArr* srcA, const CvArr* srcB,
                                      CvArr* velx, CvArr* vely );
OPENCVAPI  void  cvCalcOpticalFlowHS( const CvArr* srcA, const CvArr* srcB,
                                      int usePrevious, CvArr* velx, CvArr* vely,
OPENCVAPI  void  cvCalcOpticalFlowPyrLK( const CvArr*  imgA, const CvArr*  imgB,
                                         CvArr*  pyrA, CvArr*  pyrB,
OPENCVAPI  void  cvCalcAffineFlowPyrLK( const CvArr*  imgA, const CvArr*  imgB,
                                        CvArr*  pyrA, CvArr*  pyrB,
OPENCVAPI  void    cvUpdateMotionHistory( const CvArr* silhouette, CvArr* mhi,
OPENCVAPI  void    cvCalcMotionGradient( const CvArr* mhi, CvArr* mask, CvArr* orientation,
OPENCVAPI  double  cvCalcGlobalOrientation( const CvArr* orientation, const CvArr* mask,
                                            const CvArr* mhi, double curr_mhi_timestamp,
OPENCVAPI  CvSeq*  cvSegmentMotion( const CvArr* mhi, CvArr* seg_mask,
OPENCVAPI  void  cvAcc( const CvArr* image, CvArr* sum,
                        const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvSquareAcc( const CvArr* image, CvArr* sqSum,
                              const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvMultiplyAcc( const CvArr* imgA, const CvArr* imgB, CvArr* acc,
                                const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI  void  cvRunningAvg( const CvArr* imgY, CvArr* imgU, double alpha,
                               const CvArr* mask CV_DEFAULT(NULL) );
OPENCVAPI int  cvCamShift( const CvArr* imgProb, CvRect  windowIn,
OPENCVAPI int  cvMeanShift( const CvArr* imgProb, CvRect  windowIn,
OPENCVAPI CvSeq* cvConvexHull2( const CvArr* input,
OPENCVAPI  int  cvCheckContourConvexity( const CvArr* contour );
OPENCVAPI CvSeq*  cvConvexityDefects( const CvArr* contour, const CvArr* convexhull,
OPENCVAPI CvBox2D cvFitEllipse2( const CvArr* points );
OPENCVAPI  void  cvCalcArrHist( CvArr** arr, CvHistogram* hist,
                                const CvArr* mask CV_DEFAULT(NULL) );
                             const CvArr* mask CV_DEFAULT(NULL) );
                             int doNotClear, const CvArr* mask )
    cvCalcArrHist( (CvArr**)img, hist, doNotClear, mask );
OPENCVAPI  void  cvCalcArrBackProject( CvArr** img, CvArr* dst,
#define  cvCalcBackProject(img, dst, hist) cvCalcArrBackProject((CvArr**)img, dst, hist)
OPENCVAPI  void  cvCalcArrBackProjectPatch( CvArr** img, CvArr* dst, CvSize range,
     cvCalcArrBackProjectPatch( (CvArr**)img, dst, range, hist, method, normFactor )
OPENCVAPI  void  cvDistTransform( const CvArr* src, CvArr* dst,
OPENCVAPI  void  cvThreshold( const CvArr*  src, CvArr*  dst,
OPENCVAPI  void  cvAdaptiveThreshold( const CvArr* src, CvArr* dst, double maxValue,
OPENCVAPI  void  cvFloodFill( CvArr* arr, CvPoint seedPoint,
                              CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI  void  cvCanny( const CvArr* src, CvArr* dst, double low_threshold,
OPENCVAPI void cvPreCornerDetect( const CvArr* src, CvArr* dst,
OPENCVAPI void  cvCornerEigenValsAndVecs( const CvArr* src, CvArr* eigenvv,
OPENCVAPI void  cvCornerMinEigenVal( const CvArr* src, CvArr* eigenval,
OPENCVAPI  void  cvFindCornerSubPix( const CvArr* src,CvPoint2D32f*  corners,
OPENCVAPI void  cvGoodFeaturesToTrack( const CvArr* image, CvArr* eig_image,
                                       CvArr* temp_image, CvPoint2D32f* corners,
                                       const CvArr* mask CV_DEFAULT(NULL));
OPENCVAPI  CvSeq*  cvHoughLines2( CvArr* image, void* line_storage, int method, 
OPENCVAPI  void  cvFitLine( const CvArr* points, CvDisType dist, double param,
OPENCVAPI  void  cvImgToObs_DCT( const CvArr* arr, float* obs, CvSize dctSize,
OPENCVAPI  void  cvKMeans2( const CvArr* samples, int cluster_count,
                            CvArr* cluster_idx, CvTermCriteria termcrit );
OPENCVAPI  void  cvUnDistortOnce( const CvArr* srcImage, CvArr* dstImage,
OPENCVAPI  void  cvUnDistortInit( const CvArr* srcImage, CvArr* undistMap,
OPENCVAPI  void  cvUnDistort( const CvArr* srcImage, CvArr* dstImage,
                              const CvArr* undistMap, int interpolate CV_DEFAULT(1));
OPENCVAPI  void  cvConvertMap( const CvArr* srcImage, const CvArr* flUndistMap,
                               CvArr* undistMap, int iterpolate CV_DEFAULT(1) );
OPENCVAPI  int  cvFindChessBoardCornerGuesses( const CvArr* arr, CvArr* thresh,

};
