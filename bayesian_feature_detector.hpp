#ifndef BAYESIAN_FEATURE_DETECTOR_HPP
#define BAYESIAN_FEATURE_DETECTOR_HPP

#include <vector>
#include <values.h> // MAX_DOUBLE
#include <cv/matrix.hh>
#include <cv/gaussian_nd.hpp>
#include <cv/image.h>
#include <cv/cvext.h>

class BayesianFeatureDetector {
	int m_nscales;
	int m_width;
	std::vector< GaussianND * > m_pos;
	std::vector< GaussianND * > m_neg;
    std::vector< std::vector< CvMat *> > m_pos_data;
    std::vector< std::vector< CvMat *> > m_neg_data;
	IplImage ** m_pyr;
	cv::Image<float> m_probeImage;

public:
	BayesianFeatureDetector( int nscales, int object_width, CvSize input_size):
		m_nscales( nscales ),
		m_width( object_width )
	{
		for(int i=0; i<nscales; i++){
			m_pos.push_back( new GaussianND( object_width*object_width ) );
			m_neg.push_back( new GaussianND( object_width*object_width ) );
			m_pos_data.push_back( std::vector<CvMat *>() );
			m_neg_data.push_back( std::vector<CvMat *>() );
		}
		m_pyr = cvCreateImagePyr( input_size, 8, 1, nscales );
	}
	BayesianFeatureDetector( const char * fname ){
		assert( this->load( fname ) );
	}

	void add_example( IplImage * im, CvPoint2D32f p ){
		CvMat * imagef;
		CvPoint2D32f scaled_pos = p;
		cvPyrDownAll(im, (CvArr **) m_pyr);
        for(int j=0; j<m_nscales; j++){
            imagef = cvCreateMat(m_width, m_width, CV_32F);
            // positive examples
            cvRotateZoomImage( m_pyr[j], imagef, scaled_pos, 0, 1 );
            m_pos_data[j].push_back(imagef);
            //cvImageSC( "pos",  imagef );
            //cvImageSC( "pos_mean", pos_mean[j] );

            // negative example -- offset from eye position by 5+[0-10] px
            for(int k=0; k<5; k++){
                CvPoint2D32f neg_pt = scaled_pos;
                float ang = drand48()*2*M_PI;
                float rad = drand48()*20+10;
                neg_pt.x += rad*cos(ang);
                neg_pt.y += rad*sin(ang);

                imagef = cvCreateMat(m_width, m_width, CV_32F);

                cvRotateZoomImage( m_pyr[j], imagef, neg_pt, 0, 1 );

                m_neg_data[j].push_back(imagef);
            }
            //cvImageSC( "neg",  imagef );
            //cvImageSC( "neg_mean", neg_mean[j] );
            //cvWaitKey(33);

            // scale eye location to next pyramid level
            if(m_pyr[j+1]){
                scaled_pos.x = scaled_pos.x * m_pyr[j+1]->width/m_pyr[j]->width;
                scaled_pos.y = scaled_pos.y * m_pyr[j+1]->height/m_pyr[j]->height;
            }
        }
	}
	
	/** calculate mean and covariance matrices for examples given */
	void train(float lambda){
        for(int i=0; i<m_nscales; i++){
			CvMat pos_mean, neg_mean;
			cvReshape(&m_pos[i]->mu(), &pos_mean, 1, m_width);
			cvReshape(&m_neg[i]->mu(), &neg_mean, 1, m_width);
            cvCalcCovarMatrix( (const CvArr **) &(m_pos_data[i][0]), m_pos_data[i].size(), 
					&m_pos[i]->sigma(), &pos_mean, CV_COVAR_NORMAL | CV_COVAR_SCALE);
            cvCalcCovarMatrix( (const CvArr **) &(m_neg_data[i][0]), m_neg_data[i].size(), 
					&m_neg[i]->sigma(), &neg_mean, CV_COVAR_NORMAL | CV_COVAR_SCALE);

            add_noise_covar( &m_pos[i]->sigma(), lambda );
            add_noise_covar( &m_neg[i]->sigma(), lambda );

			m_pos[i]->recalc();
			m_neg[i]->recalc();
        }		
	}

	/**
	 * dirname/scale%02d/{pos_mean, pos_covar, pos_inv_covar}
	 */
	bool load( const char * dirname ){
		CvSize sz;
		char text[256];

		snprintf(text, 256, "%s/params.txt", dirname);
		FILE * f = fopen(text, "r");
		assert(f);
		assert(fscanf(f, "%d %d %d %d", &m_nscales, &m_width, &sz.width, &sz.height)==4);
		fclose(f);
		
		for(int i=0; i<m_nscales; i++){
			snprintf(text, 256, "%s/scale%02d/pos", dirname, i);
			m_pos.push_back( new GaussianND(text) );
			assert(m_pos.back());
			snprintf(text, 256, "%s/scale%02d/neg", dirname, i);
			m_neg.push_back(  new GaussianND(text) );
			assert(m_neg.back());
		}

		m_pyr = cvCreateImagePyr( sz, 8, 1, m_nscales );

		return m_nscales > 0;
	}

	/**
	 * dirname/scale%02d/{pos_mean, pos_covar, pos_inv_covar}
	 */
	bool save( const char * fname ){
		char text[256];
		if(mkdir( fname, 0777 )!=0 && errno!=EEXIST){
			perror("mkdir");
			return false;
		}

		snprintf(text, 256, "%s/params.txt", fname);
		FILE * f = fopen(text, "w");
		assert(f);
		fprintf(f, "%d %d %d %d\n", m_nscales, m_width, m_pyr[0]->width, m_pyr[0]->height);
		fclose(f);

		for(int i=0; i<m_nscales; i++){
			snprintf(text, 256, "%s/scale%02d", fname, i);
			if(mkdir( text, 0777 )!=0 && errno!=EEXIST){
				perror("mkdir");
				return false;
			}
			snprintf(text, 256, "%s/scale%02d/pos", fname, i);
			m_pos[i]->save(text);
			snprintf(text, 256, "%s/scale%02d/neg", fname, i);
			m_neg[i]->save(text);
		}
		return true;
	}

	/** clear examples given so far */
	void clear(){
		for(int i=0; i<m_nscales; i++){
			for(size_t j=0; j<m_pos_data[i].size(); j++){
				cvReleaseMat( &(m_pos_data[i][j]) );
			}
			for(size_t j=0; j<m_neg_data[i].size(); j++){
				cvReleaseMat( &(m_neg_data[i][j]) );
			}
			
			m_pos_data[i].clear();
			m_neg_data[i].clear();
		}
	}

	/** find point in image that maximizes log-likelihood ratio */
	CvPoint2D32f find( IplImage * im ){
		m_probeImage.realloc( *m_pyr[m_nscales-1] );
		cvPyrDownAll( im, (CvArr **) m_pyr );

		// start with coarsest level
		CvPoint origin_pt = cvPoint(0,0);
		CvPoint pt;
		for(int j=m_nscales-1; j>=0; j--){

			// fix origin to not go over image boundaries
			origin_pt.x = MAX(0, MIN(origin_pt.x, m_pyr[j]->width-m_probeImage.width));
			origin_pt.y = MAX(0, MIN(origin_pt.y, m_pyr[j]->height-m_probeImage.height));

			cvSetImageROI( m_pyr[j], cvRect( origin_pt.x, origin_pt.y, m_probeImage.width, m_probeImage.height) );
			cvScale( m_pyr[j], &m_probeImage );
			cvResetImageROI( m_pyr[j] );

			find_match( &m_probeImage, &m_pos[j]->mu(), &m_pos[j]->sigma_inv(), 
					&m_neg[j]->mu(), &m_neg[j]->sigma_inv(), &pt );

			pt.x += origin_pt.x;
			pt.y += origin_pt.y;

			if(j>0){
				origin_pt.x = pt.x*m_pyr[j-1]->width/m_pyr[j]->width;
				origin_pt.y = pt.y*m_pyr[j-1]->height/m_pyr[j]->height;
			}
		}
		return cvPoint2D32f( pt.x+0.5*m_width, pt.y+0.5*m_width );
	}	
	void add_noise_covar( CvArr * covar, float lambda ){
		CvMat diag;
		cvGetDiag( covar, &diag );
		cvAddS( &diag, cvScalar(lambda), &diag );
	}
	double gauss_log_pdf( CvArr * x, CvArr * mean, CvArr * covar_inv, CvArr * tmp1, CvArr * tmp2){
		CvMat tmp1_stub;
		CvMat tmp2_stub;
		cvSub(x, mean, tmp1);

		tmp1 = cvReshape( tmp1, &tmp1_stub, 1, 1);
		tmp2 = cvReshape( tmp2, &tmp2_stub, 1, 1);
		cvMatMul( tmp1, covar_inv, tmp2 );
		return cvDotProduct( tmp1, tmp2 );
	}
	void find_match( IplImage * data, CvArr * meanarr, CvMat * covar_inv, 
			         CvArr * neg_meanarr, CvMat * neg_covar_inv, CvPoint * res){
		CvMat mean_stub, *mean = cvGetMat(meanarr, &mean_stub);
		CvMat neg_mean_stub, *neg_mean = cvGetMat(neg_meanarr, &neg_mean_stub);

		cvReshape( meanarr, mean, 1, cvFloor(sqrt(mean->rows*mean->cols)) );
		cvReshape( neg_meanarr, neg_mean, 1, mean->rows );

		// slide patch along data and calculate (x-mean)*inv_covar*(x-mean)
		int xmax = data->width-mean->cols;
		int ymax = data->height-mean->rows;
		double err, min_err = MAXDOUBLE;
		CvMat * tmp1 = cvCreateMat(mean->rows, mean->cols, mean->type);
		CvMat * tmp2 = cvCreateMat(cvGetNumElements(mean), 1, CV_MAT_DEPTH(mean->type));

		cv::Image<float> resImage(xmax, ymax);

		//resImage = *data;
		for(int j=0; j<ymax; j++){
			for(int i=0; i<xmax; i++){
				cvSetImageROI( data, cvRect(i,j,mean->cols,mean->rows) );
				err = gauss_log_pdf( data, mean, covar_inv, tmp1, tmp2 )
				//	- log(cvGetReal2D(prior, j+y2, i+x2))
					- gauss_log_pdf( data, neg_mean, neg_covar_inv, tmp1, tmp2 ); 
				// +log(cvGetReal2D(neg_prior, j+y2, i+x2));
				//err *= log(cvGetReal2D( prior, j+y2, i+x2 ));
				//resImage[j][i] = err;
				//err -= log(cvGetReal2D(prior, j, i)/cvGetReal2D(neg_prior, j, i));
				//printf("%f ", resImage[j][i]);
				cvResetImageROI( data );
				if(err<min_err){
					//printf("%d %d\n", i, j);
					res->x = i;
					res->y = j;
					min_err = err;
					//cvCircle( &resImage, cvPoint(i+mean->width/2, j+mean->height/2), 1, CV_RGB(255,255,255));
					//cvRectangle( &resImage, cvPoint(i,j), cvPoint(i+mean->width, j+mean->height), CV_RGB(255,   255,255), 1);
				}
				//resImage[j][i] = cvGetReal2D(prior, j, i);//err;
			}
			//printf("\n");
		}
		//resImage.imagesc("win");
		//cvWaitKey();
		cvReleaseMat( &tmp1 );
		cvReleaseMat( &tmp2 );
	}
};

#endif // BAYESIAN_FEATURE_DETECTOR_HPP

