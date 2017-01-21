#include <cxcore.h>
#include <cv.h>
#include <cv/cvext.h>
#include <cv/cvransac.h>

IplImage * im;

void cvAAMDrawShape( CvArr * im, CvMat * shape, CvScalar color, bool nonnegative ) {
	CvMat m;
	CvSize sz = cvGetSize( im );
	if(shape->rows == 1){
		// interpret as Mx1 2ch matrix
		shape = cvReshape(shape, &m, 2, shape->rows*shape->cols*CV_MAT_CN(shape->type)/2);
	}
	float s = 2;
	for(int i=0; i<shape->rows; i++){
		CvPoint2D32f p = * ((CvPoint2D32f *) cvPtr1D(shape, i));
		if(nonnegative){
			p.x *= sz.width;
			p.y *= sz.height;
		}
		else{
			p.x = p.x*sz.width/2*s + sz.width/2;
			p.y = p.y*sz.height/2*s + sz.height/2;
		}
        cvCircle(im, cvPointFrom32f( p ), 2, color);
    }
}

void icvDrawPoint( CvArr * im, CvPoint2D32f p, CvScalar color ){
	CvSize sz = cvGetSize( im );
	p.x = p.x*sz.width/2 + sz.width/2;
	p.y = p.y*sz.height/2 + sz.height/2;
	cvCircle(im, cvPointFrom32f( p ), 2, color );
}

// here shape, refShape and outShape are a sequences of CV_32F points 
// can be dimension npoints x 2, 2ch x npoints x 1, etc
double icvCalcShapeRotation(CvMat * shape, CvMat * refShape){
    CvMat mRef, mS;

    // interpret data as matrices ( nbPoints x 2 columns )    
    cvReshape( shape, &mS, 1, shape->cols*shape->rows*CV_MAT_CN(shape->type)/2 );
    cvReshape( refShape, &mRef, 1, shape->cols*shape->rows*CV_MAT_CN(shape->type)/2 );

    // calculate the rotation by minimizing the sum of squared
    // point distances
    CvMat * res = cvCreateMat( 2, 2, CV_MAT_TYPE( mRef.type ) );
    CvMat * S = cvCreateMat( 2, 2, CV_MAT_TYPE( mRef.type ) );
    CvMat * V = cvCreateMat( 2, 2, CV_MAT_TYPE( mRef.type ) );
    CvMat * U = cvCreateMat( 2, 2, CV_MAT_TYPE( mRef.type ) );

    // res = mRef^T * mS
    cvGEMM( &mRef, &mS, 1, NULL, 0, res, CV_GEMM_A_T );

    cvSVD( res, S, U, V );

    // res = V*U.Transposed();
    cvGEMM( V, U, 1, NULL, 0, res, CV_GEMM_B_T );

    // res holds now a normal 2x2 rotation matrix
    double angle;
    double cos_theta = cvGetReal2D(res, 0, 0);
    double sin_theta = cvGetReal2D(res, 1, 0);
    const double epsilon = 1e-12;

    if ( 1.0-cos_theta<epsilon ) {
        // cos_theta=1  =>  shapes are already aligned
        angle = 0;
        return angle;
    }
    if ( fabs(cos_theta)<epsilon ) {
        // cos_theta=0  =>  90 degrees rotation
        return M_PI/2;
    }

    if ( 1.0+cos_theta<epsilon ) {
        // cos_theta=-1  =>  180 degrees rotation
        angle = M_PI;
    } 
	else {
        // get the rotation in radians
        double a_cos = acos( cos_theta );
        double a_sin = asin( sin_theta );

        if (a_sin<0) {
            // lower half of the unit circle
            angle = -a_cos;
        } else {
            // upper half of the unit circle
            angle = a_cos;
        }
    }

    cvReleaseMat( &res );
    cvReleaseMat( &S );
    cvReleaseMat( &U );
    cvReleaseMat( &V );
    return angle;
}

// here shape, refShape and outShape are a sequences of CV_32F points 
// can be dimension npoints x 2, 2ch x npoints x 1, etc
void icvRotateShape(CvMat * shape, double theta, CvMat * outShape){

    // interpret as 2d matrix
    CvMat shape2, outshape2;
    // interpret data as matrices ( nbPoints x 2 columns )    
    cvReshape( shape, &shape2, 1, shape->cols*shape->rows*CV_MAT_CN(shape->type)/2 );
    cvReshape( outShape, &outshape2, 1, shape->cols*shape->rows*CV_MAT_CN(shape->type)/2 );

    // set up rotation matrix
    double c00 =  cos( theta );
    double c01 = -sin( theta );
    double c10 =  sin( theta );
    double c11 =  cos( theta );

    CvMat * R = cvCreateMat( 2,2, CV_MAT_TYPE( shape2.type ) );
    cvSetReal2D( R, 0, 0, c00);
    cvSetReal2D( R, 0, 1, c01);
    cvSetReal2D( R, 1, 0, c10);
    cvSetReal2D( R, 1, 1, c11);


    cvGEMM(&shape2, R, 1, NULL, 0, &outshape2, CV_GEMM_B_T);

    cvReleaseMat(&R);

}

// here shape, refShape and outShape are a sequences of CV_32F points 
// can be dimension npoints x 2, 2ch x npoints x 1, etc
void cvShapeAlignTo( CvMat * shape, CvMat * refShape, CvMat * outShape, float * pTheta){
    // assume shape and refShape are already normalized

    // align rotation between this and refCpy   
    double theta;
    theta = icvCalcShapeRotation( shape, refShape );
    icvRotateShape( shape, theta, outShape );

    if (pTheta) {
        *pTheta = -theta;
    }
}

int icvAlignFitFunc (CvMat * data, CvMat * indices, CvMat * model);
int icvAlignDistFunc (CvMat * data, CvMat * model, double dist, CvMat * inliers);
int icvAlignDegenFunc (CvMat * data, CvMat * indices);

void cvNormalizeShape( CvMat * shape, CvMat * outShape, CvScalar * outCenter, double * outScale ){
	*outCenter = cvAvg( shape );
	cvSubS( shape, *outCenter, outShape );
	*outScale = cvNorm( outShape );
	cvScale( outShape, outShape, 1.0/(*outScale), 0);
}

void cvTransformShape( CvMat * shape, CvMat *outShape, float dx, float dy, float ds, float theta ){
	for(int i=0; i<shape->rows; i++){
		CvPoint2D32f p = *(CvPoint2D32f *) cvPtr1D( shape, i );
		CvPoint2D32f * q = (CvPoint2D32f *) cvPtr1D( outShape, i);
		q->x = ds * ( cos(theta) * p.x - sin(theta)*p.y ) + dx;
		q->y = ds * ( sin(theta) * p.x + cos(theta)*p.y ) + dy;
	}
}

void icvAlignGetModel(CvMat * model, float * dx, float * dy, float * ds, float * dtheta){
    *dx = *(float *) CV_MAT_ELEM_PTR( *model, 0, 0 );
    *dy = *(float *) CV_MAT_ELEM_PTR( *model, 1, 0 );
    *ds = *(float *) CV_MAT_ELEM_PTR( *model, 2, 0 );
    *dtheta = *(float *) CV_MAT_ELEM_PTR( *model, 3, 0 );
}

void icvAlignSetModel(CvMat * model, float dx, float dy, float ds, float dtheta){
	*(float *) CV_MAT_ELEM_PTR( *model, 0, 0 ) = dx;
	*(float *) CV_MAT_ELEM_PTR( *model, 1, 0 ) = dy;
	*(float *) CV_MAT_ELEM_PTR( *model, 2, 0 ) = ds;
	*(float *) CV_MAT_ELEM_PTR( *model, 3, 0 ) = dtheta;
}

// robustly calculate best matching rotation, translation, scaling and indicate outliers
void cvAlignShape( CvMat * shape, CvMat * refShape, float *dx, float *dy, float *dtheta, float *ds, CvMat * outInliers){
	double inlierdist = 0.1;
	int ntrials = 50;
	int nfit = 10;
	int maxdegentrials = 5;
	int dim = CV_MAT_CN(shape->type)*shape->rows*shape->cols/2;
	CvMat * data = cvCreateMat(dim, 2, CV_32FC2);
	CvMat * model = cvCreateMat(4, 1, CV_32F);
	CvMat * inliers = outInliers;

	CvMat col, colmat;
	// normalize, recenter shape 
	cvGetCol(data, &col, 0);
	cvReshape(shape, &colmat, 2, dim);  // access shape as a CV_32FC2 Mx1 matrix 
	//cvNormalizeShape( &colmat, &col );
	cvCopy(&colmat, &col);
	
	// normalize, recenter refShape 
	cvGetCol(data, &col, 1);
	cvReshape(refShape, &colmat, 2, dim);  // access shape as a CV_32FC2 Mx1 matrix 
	//cvNormalizeShape( &colmat, &col );
	cvCopy(&colmat, &col);

	// this is a ransac optimization
	cvRansac( data, inliers, model, icvAlignFitFunc, icvAlignDistFunc, icvAlignDegenFunc, ntrials, nfit, inlierdist, maxdegentrials );

	// output parameters
	icvAlignGetModel( model, dx, dy, ds, dtheta );
}


// each row of data is a correspondance pair 
// model is dx,dy,dtheta,ds
// given indices -- calculate model
int icvAlignFitFunc (CvMat * data, CvMat * indices, CvMat * model){
	CvMat * fitdata = cvCreateMat( indices->rows, data->cols, CV_32FC2 );
	float dx,dy,ds,dtheta;
	int idx;

	for(int i=0; i<indices->rows; i++){
		CvMat srow,drow;
		idx = CV_MAT_ELEM( *indices, int, i, 0 );
		cvGetRow( data, &srow, idx );
		cvGetRow( fitdata, &drow, i );
		cvCopy( &srow, &drow );
	}

	CvMat shape, refShape;
	CvScalar c1, refC;
	double scale1, refScale;

	cvZero(im);
	cvGetCol( fitdata, &shape, 0 );
	//cvAAMDrawShape( im, &shape, CV_RGB(255,255,0), false);
	cvNormalizeShape( &shape, &shape, &c1, &scale1);
	cvGetCol( fitdata, &refShape, 1 );
	//cvAAMDrawShape( im, &refShape, CV_RGB(0, 255,255), false);
	cvNormalizeShape( &refShape, &refShape, &refC, &refScale);

	//cvShowImage("win", im);
	//cvWaitKey(-1);
	cvZero(im);

	// now calculate rotation
	float theta = icvCalcShapeRotation(&shape, &refShape); 

	ds = refScale/scale1;
	dtheta = theta;
	dx = refC.val[0] - ds * ( c1.val[0]*cos(theta) - c1.val[1]*sin(theta) );
	dy = refC.val[1] - ds * ( c1.val[0]*sin(theta) + c1.val[1]*cos(theta) );
	icvAlignSetModel( model, dx, dy, ds, dtheta );

	printf("icvAlignFitFunc: (%f %f %f) (%f %f %f)\n", c1.val[0], c1.val[1], scale1, refC.val[0], refC.val[1], refScale);

	cvTransformShape( &shape, &shape, c1.val[0], c1.val[1], scale1, 0 );
	cvTransformShape( &refShape, &refShape, refC.val[0], refC.val[1], refScale, 0 );

	cvAAMDrawShape( im, &shape, CV_RGB(255,255,0), false);
	cvAAMDrawShape( im, &refShape, CV_RGB(0,255,255), false);
	cvTransformShape( &shape, &shape, dx, dy, ds, 0 );
	cvShowImage("win", im);
	//cvWaitKey(-1);
	cvAAMDrawShape( im, &shape, CV_RGB(255,255,0), false);
	cvShowImage("win", im);
	//cvWaitKey(-1);
	

	return 0;
}


// calculate L2 distance between point pairs and count inliers
int icvAlignDistFunc (CvMat * data, CvMat * model, double dist, CvMat * inliers){
	CvPoint2D32f p, refP, p3;
	float dx,dy,dtheta,ds;
	icvAlignGetModel( model, &dx, &dy, &ds, &dtheta );
	
	int ninliers = 0;
	cvZero(im);
	for(int i=0; i<data->rows; i++){
		p = *(CvPoint2D32f *) cvPtr2D(data, i, 0);
		refP = *(CvPoint2D32f *) cvPtr2D(data, i, 1);

		// transform p
		p3.x = p.x*cos(dtheta)*ds - p.y*sin(dtheta)*ds + dx;
		p3.y = p.x*sin(dtheta)*ds + p.y*cos(dtheta)*ds + dy;
		
		icvDrawPoint( im, refP, CV_RGB(0,0,255));
		
		if(L2(p3, refP)<dist){
			if(inliers) cvSet1D(inliers, i, cvScalarAll(1.0) );
			icvDrawPoint( im, p3, CV_RGB(0,255,0));
			ninliers++;
		}
		else{
			if(inliers) cvSet1D(inliers, i, cvScalarAll(0.0) );
			icvDrawPoint( im, p3, CV_RGB(255,0,0));
		}
	}
	cvShowImage("win", im);
	cvWaitKey(-1);
	printf("ninliers: %d\n", ninliers);
	return ninliers;
}

// always return 0 for now
int icvAlignDegenFunc (CvMat * data, CvMat * indices){
	return 0;
}

int main(int argc, char ** argv){
	CvMat * pts = cvCreateMat( 20, 1, CV_32FC2 );
	CvMat * cpts = cvCreateMat( 20, 1, CV_32FC2 );
	cvNamedWindow("win", 0);
	im = cvCreateImage(cvSize(250,250),8,3);
	cvRedirectError( cvSegReport );

	// create a set of points in a circle
	for(int i=0; i<pts->rows; i++){
		float theta = M_PI*i/pts->rows;
		CvPoint2D32f * p = (CvPoint2D32f *) CV_MAT_ELEM_PTR( *pts, i, 0);
		p->x = cos(theta)*.25;
		p->y = sin(theta)*.25;
	}

	// now scale, translate and corrupt points
	float s = .8 + .1 * drand48();
	float theta = M_PI*drand48();
	CvPoint2D32f trans = cvPoint2D32f( drand48(), drand48() );
	for(int i=0; i<pts->rows; i++){
		CvPoint2D32f p = CV_MAT_ELEM( *pts, CvPoint2D32f, i, 0);
		CvPoint2D32f * q = (CvPoint2D32f *) CV_MAT_ELEM_PTR(*cpts, i, 0 );
		if(drand48()>0.25){
			q->x = cos(theta)*s*p.x - sin(theta)*s*p.y + trans.x + (drand48()-0.5)*.1;
			q->y = sin(theta)*s*p.x + cos(theta)*s*p.y + trans.y + (drand48()-0.5)*.1;
		}
		else{
			q->x = .5*(drand48()-0.5);
			q->y = .5*(drand48()-0.5);
		}
		
	}

	// draw them
	float dx,dy,dtheta,ds;
	cvAlignShape( cpts, pts, &dx, &dy, &dtheta, &ds, NULL);

	cvZero(im);

	cvAAMDrawShape( im, pts, CV_RGB(0, 255, 0) , false);
	cvAAMDrawShape( im, cpts, CV_RGB(255, 0, 0), false);
	
	cvTransformShape( cpts, cpts, dx, dy, ds, 0);
	cvAAMDrawShape( im, cpts, CV_RGB(0, 0, 255), false);
	printf("Solved: %f %f %f %f\n", dx, dy, dtheta, ds);
	printf("Actual: %f %f %f %f\n", -(trans.x*cos(theta)+trans.y*sin(theta))/s, (trans.x*sin(theta)-trans.y*cos(theta))/s, -theta, 1/s);
	cvShowImage("win", im);
	cvWaitKey(-1);
	
}
