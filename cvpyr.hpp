#ifndef __CV_IMAGE_PYRAMID_HH
#define __CV_IMAGE_PYRAMID_HH

struct CvImagePyr {
	IplImage ** images;
	int size;
};

CvImagePyr * cvCreateImagePyr( int nimages, CvSize size, int depth, int nch, int minSize=16, int maxLevels=0 ){
	if(maxLevels=0){
        // calculate number of scales to consider
		int logMinSize = log2(minSize);
        int w = (int) log2(size.width);
        int h = (int) log2(size.height);
        if(w>h) maxLevels = h - logMinSize + 1;
        else maxLevels = w - logMinSize + 1;
	}

	CvImagePyr * pyr = cvAlloc( sizeof( CvImagePyr ) );
	pyr->images = cvAlloc (sizeof((IplImage *) * maxLevels));
	pyr->size = maxLevels;
	
	for(int i=0; i<maxLevels; i++){
		pyr->images[i] = cvCreateImage(size, depth, nch);
		size.width /= 2;
		size.height /= 2;
	}

	return pyr;
}

void cvPyrDownAll( CvImagePyr * pyr ) {
	for(int i=1; i<pyr->size; i++){
		cvPyrDown(pyr->images[i-1], pyr->images[i], CV_GAUSSIAN_5x5 );
	}
}

void cvPyrUpAll( CvImagePyr * pyr ) {
	for(int i=pyr->size-1; i>0; i++){
		cvPyrUp(pyr->images[i], pyr->images[i-1], CV_GAUSSIAN_5x5 );
	}
}

#endif //__CV_IMAGE_PYRAMID_HH
