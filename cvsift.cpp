#include <cxcore.h>
#include <cv.h>
#include <cv/cvext.h>
#include <matrix.hh>
#include "cvsift.h"

#define CV_SIFT_RBIN( nbins, angle ) ( cvRound(nbins*angle*.5/M_PI) )

CvSIFTDescriptor * cvCreateSIFTDescriptor(int sbins, int rbins){
	CvSIFTDescriptor * s;
	s = (CvSIFTDescriptor *) cvAlloc( sizeof(CvSIFTDescriptor) );
	cvInitMatHeader( &s->hist, 1, sbins*sbins*rbins, CV_32F );
	cvCreateData( &s->hist);
	s->sbins = sbins;
	s->rbins = rbins;
	return s;
}

/** Allocate scale space data structure */
CvScaleSpacePyr * cvCreateScaleSpacePyr( CvSize size, int nscales, int minSize, int maxLevels ){
	CvScaleSpacePyr * pyr;
	
	if(nscales<1) {
		nscales = 1;
	}

	// for n sub scales, need n+3 gaussian images
	// 1 extra to compute DoG 
	// 2 extra for searching for min/max

        // for k levels in the scale space,
        // the freqFactor between them is 2^1/k
        // note that to search all scales, we need overlap
        // between scale space octaves.

        // Normally, an octave would have scales sigma -> 2^(k-1)/k*sigma 
        // However, in this case, not all scales are searched.
        //
        // scale       gauss  dog   search
        // s            Y      N       N
        // s*2^1/k      Y      Y       N
        // s*2^2/k      Y      Y       Y
        // ...          Y      Y       Y
        // s*2^(k-1)/k  Y      Y       N

        // Specifically, the first two scales don't get searched, and
        // the last doesn't.
        // We make up for this by adding some overlap between octaves

        // ...
        // s*2^(k-1)/k  Y      Y       Y
        // s*2          Y      Y       Y   // also in octave 1
        // s*2*2^1/k    Y      Y       Y   // also in octave 1
        // s*2*2^1/k    Y      Y           // also in octave 1

        // Now s*2^n and s*2^(n+1/k) are searched (except for 2^0 and 2^1/k)
        // and s*2^(n+(k-1)/k) is searched (except for the very end of the pyramid)

	pyr = (CvScaleSpacePyr *) cvAlloc( sizeof(CvScaleSpacePyr) );

	// calculate # of levels
	pyr->nLevels = MIN( (log( MIN(size.width, size.height) ) - log(minSize) )/log( 2.0 ), maxLevels ); 
	printf("Size %dx%d nLevels %d\n", size.width, size.height, pyr->nLevels );

	 pyr->levels = (CvScaleSpace *) cvAlloc( sizeof(CvScaleSpace) * pyr->nLevels );
	for(int i=0; i<pyr->nLevels; i++){
		assert(size.height>=minSize);
		assert(size.width>=minSize);
		pyr->levels[i].nScales = nscales;
		pyr->levels[i].images = (CvMat **)cvAlloc( sizeof(void*) * nscales);
		for(int j=0; j<nscales; j++){
			pyr->levels[i].images[j] = cvCreateMat( size.height, size.width, CV_32F );
		}
		size.height/=2;
		size.width/=2;
	}
	return pyr;
}

void cvReleaseScaleSpacePyr( CvScaleSpacePyr ** scalespace_ptr ){
	CvScaleSpacePyr * scalespace = *scalespace_ptr;
	for(int i=0; i<scalespace->nLevels; i++){
		for(int j=0; j<scalespace->levels[i].nScales; j++){
			cvReleaseMat( &scalespace->levels[i].images[j] );
		}
		cvFree( &(scalespace->levels[i].images ) );
	}
	cvFree( &(scalespace->levels) );
	cvFree( scalespace_ptr );
}

/** Compute SIFT descriptors for scale space features */
CvSeq * cvComputeSIFT( CvScaleSpacePyr * scalespace, CvSeq * features, CvMemStorage storage ){

}

CvSIFTDescriptor * cvGetSeqElemSIFT( CvSeq * seq, int idx ){

}


void cvCalcGaussianOctave( const CvArr * src, CvArr ** dst, int len, float sigma, int nscales){
	float s;
	int i;
	if(len<=0){
		return;
	}
	float kTerm = sqrt( pow( exp2( 1.0/nscales), 2.0 ) - 1.0 );

	for(i=0; i<len; i++){
		s = sigma * exp2(i/(double)nscales) * kTerm; // s * 2^i/len
		printf("sigma = %f\n", s);
		cvSmooth( src, dst[i], CV_GAUSSIAN, 0, 0, s );
		//src=dst[i]; // c# implementation uses previously smoothed image ???
	}
}

/** Compute all levels of scale space */
void cvComputeScaleSpacePyr( const CvArr * image, CvScaleSpacePyr * scalespace, float sigma ){
	int nscales = scalespace->levels[0].nScales;
	/*			for(int y=0; y<3; y++){
				for(int x=0; x<3; x++){
					printf("%f ", cvGetReal2D(image,y,x) );
				}
				printf("\n");
			}
	*/

	cvResize(image, scalespace->levels[0].images[0], CV_INTER_CUBIC);
	// (Lowe03, p10, "We assume that the original image has a blur of at
	// least \sigma = 0.5 ...")
	// Then from iLab Neuromorphic Vision C++ Toolkit
	// saliency/src/SIFT/VisualObject.C
	// To feed the first ScaleSpace in our series, apply some
	// initial blur so that the input image has an effective blur of
	// the desired sigma. We assume that the original image has a
	// blur of at least 0.5 by construction. Since its size has been
	// doubled (octscale=0.5), then that becomes 1.0. We assume that
	// the sigma=1.6 applies to the doubled image. Remember that the
	// variances add when we sequentially convolve by
	// Gaussians. Hence the additional blur we need is such that
	// sigma^2 = 1^2 + blursig^2:

	float sigma0 = sqrt(sigma * sigma - 1.0F);
	cvSmooth( scalespace->levels[0].images[0], scalespace->levels[0].images[0], CV_GAUSSIAN, 0, 0, sigma0 );

	for(int i=0; i<scalespace->nLevels; i++){
		if(i>0){
			// second to last gaussian map represents half-scale
			cvResize( scalespace->levels[i-1].images[nscales-3], 
					  scalespace->levels[i].images[0], CV_INTER_NN );
		}
		cvCalcGaussianOctave( scalespace->levels[i].images[0], (CvArr **) scalespace->levels[i].images+1, 
				scalespace->levels[i].nScales-1, sigma, nscales-3);
		/*for(int j=0; j<nscales; j++){
			printf("Image scale %d\n", j);
			for(int y=0; y<3; y++){
				for(int x=0; x<3; x++){
					printf("%f ", cvGetReal2D(scalespace->levels[i].images[j],y,x) );
				}
				printf("\n");
			}
		}*/
	}
}

#if 0
void cvCalcDogOctave( const CvArr ** gauss_octave, CvArr ** dog_octave, int len){
	int i=0;
	for(i=0; i<len; i++){
		cvSub(gauss_octave[i+1], gauss_octave[i], dog_octave[i]);
	}
}
#endif 

void cvComputeDoGScaleSpacePyr( CvScaleSpacePyr * scalespace, CvScaleSpacePyr * dogspace ){
	for(int i=0; i<scalespace->nLevels; i++){
		for(int j=1; j<scalespace->levels[i].nScales; j++){
			cvSub( scalespace->levels[i].images[j], scalespace->levels[i].images[j-1],
					dogspace->levels[i].images[j-1] );
		}
	}
}

template <typename T>
void cvLocalPeaks3( const Matrix<T> & im1, const Matrix<T> & im2, const Matrix<T> & im3, float min_val, 
						 CvSeq * seq ){   
	for(int j=1; j<im1.rows-1; ++j){
		for(int i=1; i<im1.cols-1; ++i){
			bool ismax=true;
			bool ismin=true;
			if(fabs(im1[j][i]) < min_val ) continue;
			for(int jj=-1; jj<=1; ++jj){
				for(int ii=-1; ii<=1; ++ii){
					if(ii==0 && jj==0){
						if(ismax) ismax = ismax && im1[j][i] >= im2[j][i] && 
						                       im1[j][i] >= im3[j][i]; 
						if(ismin) ismin = ismin && im1[j][i] <= im2[j][i] && 
						                       im1[j][i] <= im3[j][i];
					}
					else{
						if(ismax) ismax = ismax && im1[j][i] > im2[j+jj][i+ii] && 
						                       im1[j][i] > im3[j+jj][i+ii] &&
											   im1[j][i] > im1[j+jj][i+ii];
						if(ismin) ismin = ismin && im1[j][i] < im2[j+jj][i+ii] && 
						                       im1[j][i] < im3[j+jj][i+ii] &&
											   im1[j][i] < im1[j+jj][i+ii];
					}
				}
			}
			if(ismin || ismax){
				CvPoint p = cvPoint(i,j);
				cvSeqPush( seq, &p );
			}
		}
	}
}
void cvLocalPeaks3( const CvArr * m1, const CvArr * m2, const CvArr * m3, float min_val, CvSeq * seq ){
	cvLocalPeaks3( *(Matrix<float> *)m1, *(Matrix<float> *)m2, *(Matrix<float> *)m3, min_val, seq );
}

#include <highgui.h>
void icvDoGFeatureDetectOctave( const CvArr ** gauss_octave, CvArr ** dog_octave, int len, int level, float min_val, CvSeq * seq ){
	CvSize sz=cvGetSize(gauss_octave[0]);
	int r=sz.height; int c=sz.width; int t=CV_32F;
	CvMat * min1 = cvCreateMat(r,c,t);
	CvMat * min2 = cvCreateMat(r,c,t);
	CvMat * *max2, *max3, *tmp;

	CvMemStorage * storage = cvCreateMemStorage(0);
	CvSeq * peakseq = cvCreateSeq( CV_32SC2, sizeof(CvSeq), sizeof(CvPoint), storage );
	CvSeqReader reader;
	CvSeqWriter writer;
	
	cvStartAppendToSeq( seq, &writer );

	float scale = 0;
	tmp = cvCreateMat( sz.height, sz.width, CV_8U );
	char str[256];
	for(int k=2; k<len; ++k){ 
		printf("peak search at level %d imsize %d %d\n", k-1, sz.width, sz.height);

		cvConvert(gauss_octave[k-1], tmp);
		sprintf(str, "image_%d_%d.png", level, k-1);
		cvSaveImage(str, tmp);
		
		cvLocalPeaks3( dog_octave[k-1], dog_octave[k-2], dog_octave[k], min_val, peakseq);

		// convert cvPoints to scale points
		cvStartReadSeq( peakseq, &reader );
		for(int i=0; i<peakseq->total; ++i){
			CvPoint * p = (CvPoint *) reader.ptr;
			CvScalePoint sp = cvScalePoint( p->x, p->y, level, k-1, 0 );
			CV_WRITE_SEQ_ELEM(sp, writer);
			CV_NEXT_SEQ_ELEM( seq->elem_size, reader );
			printf("%d %d\n", p->x, p->y);
		}
		printf("found %d features at level %d\n", peakseq->total, k-1);
		cvClearSeq( peakseq );
	}
	cvEndWriteSeq( &writer );
}

CvSeq * cvFindFeaturesDoG(  const CvArr * im, CvMemStorage * storage, CvSeq ** features, 
		float dogThresh, float edgeRatio, float contrastThresh, float sigma,    
		int minImageSize, int nsub ){
	CvSize sz = cvGetSize(im);
	sz.width = 2*sz.width;
	sz.height = 2*sz.height;
	//
	CvScaleSpacePyr * scale_space = cvCreateScaleSpacePyr( sz, nsub+3, minImageSize, 5 );
	CvScaleSpacePyr * dog_scale_space = cvCreateScaleSpacePyr( sz, nsub+2, minImageSize, 5 );

	cvComputeScaleSpacePyr( im, scale_space );
	cvComputeDoGScaleSpacePyr( scale_space, dog_scale_space );

	if(*features==NULL){
		*features = cvCreateSeq( 0, sizeof(CvSeq), sizeof(CvScalePoint), storage );
	}

	//printf("FindPeaks: scale 0.50, testing %d levels\n", dog_scale_space->levels[0].nScales-2);
	for(int i=0; i<scale_space->nLevels; i++){
		icvDoGFeatureDetectOctave( (const CvArr **) scale_space->levels[i].images, (CvArr **) dog_scale_space->levels[i].images,
				 dog_scale_space->levels[i].nScales, i, dogThresh, *features);
		//printf("Found %d features thus far \n", (*features)->total);
	}

	//cvRefineFeatures( scale_space, *features, edgeRatio, contrastThresh );
	
	cvReleaseScaleSpacePyr( &scale_space );
	cvReleaseScaleSpacePyr( &dog_scale_space );

	return *features;
}

#if 0
	/*void calcGradient(const int &x, const int &y, const cv::Image<float> & space, 
	                  float & magnitude, float & angle){
		// compute gradient magnitude, orientation
		float dx = space[y][x+1]-space[y][x-1];
		float dy = space[y+1][x]-space[y-1][x];
		magnitude = sqrt( dx*dx + dy*dy );
		angle = atan2(dy,dx);
	}*/
	
	// Ix and Iy are pre rotated to orientation of point
	// feature is assumed to be located at the center of image patch
	template <typename image_t>
	void calcFromPatch(const ScalePoint &p, const image_t & Ix, const image_t & Iy) {
		(*this)=0.0;
		int n = this->length();
		int radius = n*size;
		int xmin = 0;//MAX(0, this->x-radius) - this->x + Ix.width/2;
		int ymin = 0;//MAX(0, this->y-radius) - this->y + Ix.height/2;
		int xmax = Ix.width;//MIN(width-1, this->x+radius) - this->x + Ix.width/2;
		int ymax = Ix.height;//MIN(height-1, this->y+radius) - this->y + Ix.height/2;
		for(int j=ymin; j<ymax; j++){
			int ybin = (j*size)/Ix.height;
			for(int i=xmin; i<xmax; i++){
				int xbin = (i*size)/Ix.width;
				
				float magnitude;
				float angle;
				ScalePoint::calcGradient(i,j,Ix,Iy,magnitude,angle);
				(*this)(xbin,ybin,ScalePoint::binGradient(angle,nbins)) += magnitude;
			}
		}
		this->normalize();
	}
	

	void calc(ScalePoint &p, const cv::DoGPyramid<float> & dog) {
		calc(p, dog.getIx(p.level, p.subLevel), dog.getIy(p.level, p.subLevel));
	}
	
	// Ix and Iy are standard sobel images from the pyramid, not preprocessed
	template <typename image_t>
	void calc(ScalePoint &p, const image_t & Ix, const image_t & Iy){
		(*this)=0.0;
		int n = this->length();
		int radius = n*size;
		int xmin = (int) (MAX(0, cvFloor(p.x)-radius));
		int ymin = (int) (MAX(0, cvFloor(p.y)-radius));
		int xmax = (int) (MIN(Ix.width, cvFloor(p.x)+radius));
		int ymax = (int) (MIN(Ix.height, cvFloor(p.y)+radius));
		int xoff = (int) (p.x-radius);
		int yoff = (int) (p.y-radius);
		for(int j=ymin; j<ymax; j++){
			int ybin = ((j-yoff)*size)/(radius*2);
			for(int i=xmin; i<xmax; i++){
				int xbin = ((i-xoff)*size)/(radius*2);

				float magnitude;
				float angle;
				ScalePoint::calcGradient(i,j,Ix,Iy,magnitude,angle);
				// need to account for angle relative to orientation
				angle-=p.angle;
				assert(xbin>=0 || fprintf(stderr, "xmin=%d ymin=%d xmax=%d ymax=%d xoff=%d yoff=%d p.x=%f p.y=%f radius=%d\n", xmin, ymin,xmax,ymax,xoff,yoff,p.x,p.y,radius)==0);
				assert(xbin<4 || fprintf(stderr, "xmin=%d ymin=%d xmax=%d ymax=%d xoff=%d yoff=%d p.x=%f p.y=%f radius=%d\n", xmin, ymin,xmax,ymax,xoff,yoff,p.x,p.y,radius)==0);
				assert(ybin>=0 || fprintf(stderr, "xmin=%d ymin=%d xmax=%d ymax=%d xoff=%d yoff=%d p.x=%f p.y=%f radius=%d\n", xmin, ymin,xmax,ymax,xoff,yoff,p.x,p.y,radius)==0);
				assert(ybin<4 || fprintf(stderr, "xmin=%d ymin=%d xmax=%d ymax=%d xoff=%d yoff=%d p.x=%f p.y=%f radius=%d\n", xmin, ymin,xmax,ymax,xoff,yoff,p.x,p.y,radius)==0);
				(*this)(xbin,ybin,ScalePoint::binGradient(angle,nbins)) += magnitude;
			}
		}
		this->normalize();
	}
};

} // namespace sift
#endif

	// Return adjustment (scale, y, x) on success,
	void icvCalcAdjustment (const Matrix<float> & below, 
			const Matrix<float> & current, 
			const Matrix<float> & above,
			int x, int y, float *dS, float *dx, float *dy, double * dp){

		/*Console.WriteLine ("GetAdjustment (point, {0}, {1}, {2}, out double dp
		  )",
		  level, x, y);*/
		*dp = 0.0;


		// Hessian
		float  H[9];
		CvMat H_mat = cvMat(3, 3, CV_32F, H);
		H[0] = below[y][x] - 2 * current[y][x] + above[y][x];
		H[1] = H[3] = 0.25 * (above[y+1][x] - above[y-1][x] -
				(below[y + 1][x] - below[y - 1][x]));
		H[2] = H[6] = 0.25 * (above[y][x + 1] - above[y][x - 1] -
				(below[y][x + 1] - below[y][x - 1]));
		H[4] = current[y - 1][x] - 2 * current[y][x] + current[y + 1][x];
		H[5] = H[7] = 0.25 * (current[y + 1][x+1] - current[y + 1][x - 1] -
				(current[y - 1][x + 1] - current[y - 1][x - 1]));
		H[8] = current[y][x - 1] - 2 * current[y][x] + current[y][x + 1];

		// derivative
		float d[3];
		d[0] = 0.5 * (above[y][x] - below[y][x]);             //dS
		d[1] = 0.5 * (current[y + 1][x] - current[y - 1][x]); //dy
		d[2] = 0.5 * (current[y][x + 1] - current[y][x - 1]); //dx

		CvMat d_mat = cvMat(3,1,CV_32F,d);

		// Solve: H*b = -d --> b = inv(H)*d
		float b[3];
		CvMat b_mat = cvMat(3,1,CV_32F,b);

		cvSolve(&H_mat, &d_mat, &b_mat); //= H.inv()*d;
		cvScale( &b_mat, &b_mat, -1 );

		*dp = cvDotProduct( &b_mat, &d_mat );

		*dS = b[0];
		*dx = b[1];
		*dy = b[2];
	}

	float icvSubLevelToScale( int sublevel ){
		return 0;
	}
	int icvScaleToSubLevel( float scale ){
		return 0;
	}

	bool icvSubPixelLocalize(CvScalePoint * sp, CvScaleSpace * space, float * contrast, int nadjustments=2){
		bool needToAdjust = true;
		int x = cvRound(sp->center.x);
		int y = cvRound(sp->center.y);
		int sublevel = sp->subLevel;
		//std::cerr<<"calcSubPixel Start: "<<p<<std::endl;
		while (needToAdjust) {

			// Reject points that lie on the border of scale space

			if (sublevel <= 0 || sublevel >= (space->nScales - 1))
				return false;

			if (x <= 0 || x >= (space->images[sublevel]->cols - 1))
				return false;
			if (y <= 0 || y >= (space->images[sublevel]->rows - 1))
				return false;

			float adjS,adjX,adjY;
			double dp;
			icvCalcAdjustment( *Matrix<float>::safe_cast(space->images[sublevel-1]), 
					           *Matrix<float>::safe_cast(space->images[sublevel]), 
							   *Matrix<float>::safe_cast(space->images[sublevel+1]), 
							   x, y, &adjS, &adjX, &adjY, &dp);

			printf("icvSubPixelLocalize %f %f %f\n", adjS, adjX, adjY);
			// Get adjustments and check if we require further adjustments due
			// to pixel level moves. If so, turn the adjustments into real
			// changes and continue the loop. Do not adjust the plane, as we
			// are usually quite low on planes in thie space and could not do
			// further adjustments from the top/bottom planes.
			if (fabs (adjX) > 0.5 || fabs (adjY) > 0.5 || fabs(adjS) > 0.5) {
				// Already adjusted the last time, give up
				if (nadjustments == 0) {
					//Console.WriteLine ("too many adjustments, returning");
					return false;
				}

				nadjustments -= 1;

				// Check that just one pixel step is needed, otherwise discard
				// the point
				double distSq = adjX * adjX + adjY * adjY;
				if (distSq > 2.0)
					return (false);

				x = (int)(x + adjX + 0.5);  
				y = (int)(y + adjY + 0.5); 
				sublevel = (int)(sublevel + adjS + 0.5);

				//point.Level = (int) (point.Level + adjS + 0.5);
			//	std::cerr<<"moved point by ("<<adjX<<","<<adjY<<": "<<adjS<<
				           //") to ("<<x<<","<<y<<": "<<sublevel<<")"<<std::endl;
				/*Console.WriteLine ("moved point by ({0},{1}: {2}) to ({3},{4}:
				  {5})",
				  adjX, adjY, adjS, point.X, point.Y, point.Level);*/
				continue;
			}

			sp->center.x = cvRound(x + adjX);
			sp->center.y = cvRound(y + adjY);
			sp->subLevel = cvRound(sublevel + adjS);

			*contrast = (*Matrix<float>::safe_cast(space->images[sublevel]))[y][x]+0.5*dp;
			
			//std::cerr<<"calcSubPixel: "<<p<<std::endl;
			
			//p.scale = _spaces[level].getSigma(
			//p.scale = exp2( (level + adjS)/_spaces.getNumSubScales()) * 1.6;
			//FIXME: csharp has p.value = space[point.X, point.Y] + 0.5 * dp;
			/* for procesesing with gnuplot
			 *
			 Console.WriteLine ("{0} {1} # POINT LEVEL {2}", point.X,
			 point.Y, basePixScale);
			 Console.WriteLine ("{0} {1} {2} # ADJ POINT LEVEL {3}",
			 adjS, adjX, adjY, basePixScale);
			 */

			// Check if we already have a keypoint within this octave for this
			// pixel position in order to avoid dupes. (Maybe we can move this
			// check earlier after any adjustment, so we catch dupes earlier).
			// If its not in there, mark it for later searches.
			//
			// FIXME: check why there does not seem to be a dupe at all
			//if (processesed[point.X, point.Y] != 0)
			//	return (true);

			//processesed[point.X, point.Y] = 1;

			// Save final sub-pixel adjustments.
			//PointLocalInformation local = new PointLocalInformation (adjS, adjX,
			//		adjY);
			//local.DValue = dp;
			//local.DValue = space[point.X, point.Y] + 0.5 * dp;
			//point.Local = local;

			//needToAdjust = false;
			return true;
		}

		return false;

	}

	bool icvFilterEdge(const Matrix<float> & space, int x, int y, double r){
		double D_xx, D_yy, D_xy;

		// first check bounds
		if(x<=0 || y<=0 || x>=space.width-1 || y>=space.height-1){
			return false;
		}
		
		// Calculate the Hessian H elements [ D_xx, D_xy ; D_xy , D_yy ]
		D_xx = space[y][x+1] + space[y][x-1] - 2.0 * space[y][x];
		D_yy = space[y + 1][x] + space[y - 1][x] - 2.0 * space[y][x];
		D_xy = 0.25 * ((space[y + 1][x + 1] - space[y - 1][x + 1]) -
				(space[y + 1][x - 1] - space[y - 1][x - 1]));

		// page 13 in Lowe's paper
		double TrHsq = D_xx + D_yy;
		TrHsq *= TrHsq;
		double DetH = D_xx * D_yy - (D_xy * D_xy);

		double r1sq = (r + 1.0);
		r1sq *= r1sq;

		//fprintf(stderr, "x=%d, y=%d, TrHsq * r = %f, DetH * r1sq = %f, TrHsq / DetH = %f, r1sq /r = %f\n",
		//		x, y, (TrHsq * r), (DetH * r1sq),
		//		(TrHsq / DetH) , (r1sq / r));
		
		// reject if ratio of curvatures is greater than threshold
		if ( fabs(TrHsq / DetH) > (r1sq / r)) {
			return true;
		}

		return false;
	}

	void icvCalcOrientation(CvScalePoint * sp, CvMat * ss_image){
		int nbins = 36;
		int hist[36];
		// TODO -- fill in sigma here
		//
		//double sigma =  * 3.0;
		double sigma = 1.0;
		double factor = 2.0*sigma*sigma;
		int radius = (int)(3.0*sigma / 2.0 + 0.5);
		int radiusSq = radius*radius;
		
		bzero(hist, sizeof(int)*36);
		
		int w = ss_image->cols;
		int h = ss_image->rows;
		int xmin = cvRound(MAX(1, sp->center.x-radius));
		int ymin = cvRound(MAX(1, sp->center.y-radius));
		int xmax = cvRound(MIN(w-1, sp->center.x+radius));
		int ymax = cvRound(MIN(h-1, sp->center.y+radius));

		float gradient_mag, gradient_ang;

		for(int i=xmin; i<xmax; i++){
			for(int j=ymin; j<ymax; j++){
				int relX = cvRound(i - sp->center.x);
				int relY = cvRound(j - sp->center.y);

				// only consider points in circle
				double d = relX*relX + relY*relY;
				if(d>radiusSq) continue;
				
				// gaussian weight
				double gauss_weight = exp(  -( d / factor)  );
				
				int bin;
				bin = CV_SIFT_RBIN( gradient_ang, nbins );
				hist[bin] += cvRound(gradient_mag * gauss_weight);
			}
		}
		
		// find maximum bin
		double maxv=hist[0];
		int maxb=0;
		for(int b=1; b<nbins; b++){
			if(hist[b]>maxv){
				maxv=hist[b];
				maxb=b;
			}
		}
		sp->angle = (M_PI*2*maxb)/nbins;
	}

/** Refine features.  Orients features, eliminates edge like features, 
 * features with poor contrast and performs sub pixel localization. */
void cvRefineFeatures( CvScaleSpacePyr * scalespace, CvSeq * features, float edgeRatio, float contrastThresh){
	int i=0;
	float contrast;

	for(int i=0; i<features->total; ){
		CvScalePoint *sp = CV_GET_SEQ_ELEM( CvScalePoint, features, i );
		int level=sp->level;
		int sub_level=sp->subLevel;

		//std::cout<<"Filtering feature: "<<*it<<std::endl;
		sp->center.x+=.49;
		sp->center.y+=.49;

		// Lowe uses an r value of 10
		// Here we are comparing (r+1)^2/r to the ratio of largest direction of curvature 
		// to smallest curvature, so larger values let in more edges
		//std::cout<<"Feature "<<i<<" - "<<*it<<" - ";
		if(icvFilterEdge( *Matrix<float>::safe_cast(scalespace->levels[ level ].images[ sub_level ]),
					cvRound(sp->center.x), cvRound(sp->center.y), edgeRatio)){

			cvSeqRemove( features, i );
			std::cout<<"Rejected -- too edge like\n";
		}
		// this step ought to spawn other features if the orientation is ambiguous
		else if(!icvSubPixelLocalize(sp, &(scalespace->levels[level]), &contrast)){
			cvSeqRemove( features, i );
			std::cout<<"Rejected -- could not sub pixel localize\n";
		}
		else if(contrast < contrastThresh){
			cvSeqRemove( features, i );
			std::cout<<"Rejected -- poor contrast"<<contrast<<"<"<<contrastThresh<<std::endl;
		}
		else{
			icvCalcOrientation(sp, Matrix<float>::safe_cast(scalespace->levels[ level ].images[ sub_level ]));
			i++;
		}
	}	
}

#define CV_SIFT_ELEM( desc, x, y, r ) (desc->hist.data.fl[x*desc->sbins*desc->rbins+y*desc->rbins+desc->rbins])

void icvComputeSIFT1(CvArr * gradientMag, CvArr * gradientAngle, CvScalePoint sp, CvSIFTDescriptor * d){
	int n;
	CvPoint2D32f p = sp.center;
	int radius;
	int xmin,ymin,xmax,ymax,xoff,yoff;
	int i,j;
	int xbin,ybin;
	int size;
	float angle,mag;
	int bin;
	int nbins = d->rbins;
	CvSize sz = cvGetSize(gradientMag);
	CvMat mag_mat_stub, *mag_mat=(CvMat *)gradientMag;
	CvMat ang_mat_stub, *ang_mat=(CvMat *)gradientAngle;

	if(!CV_IS_MAT(ang_mat)) ang_mat = cvGetMat( gradientAngle, &ang_mat_stub );
	if(!CV_IS_MAT(mag_mat)) mag_mat = cvGetMat( gradientAngle, &mag_mat_stub );

	cvZero(d);
	n = cvGetNumElements( d );
	radius = n*d->sbins;
	xmin = MAX(0, cvFloor( p.x-radius ));
	ymin = MAX(0, cvFloor( p.y-radius ));
	xmax = MIN(sz.width, cvFloor(p.x+radius));
	ymax = MIN(sz.height, cvFloor(p.y+radius));

	for( j=ymin; j<ymax; ++j ){
		ybin = ((j-yoff)*size)/(radius*2);
		float * mag_ptr = CV_MAT_ROW_PTR(*mag_mat, j, float);
		float * ang_ptr = CV_MAT_ROW_PTR(*ang_mat, j, float);
		for( i=xmin; i<xmax; ++i ){
			xbin = ((i-xoff)*size)/(radius*2);

			angle -= sp.angle;
			bin = CV_SIFT_RBIN( nbins, angle);
			CV_SIFT_ELEM(d, xbin, ybin, bin) += mag;
		}
	}
	cvScale(d, d, 1.0/cvNorm(d), 0);
}

//#endif //__SIFT_DESCRIPTOR_H
