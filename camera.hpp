#include <cxcore.h>
#include <cvaux.h>
#include <cv/matrix.hh>

namespace cv {

struct Camera : CvCamera {
	CvSize2D32f & size() {
		return *(CvSize2D32f *)(this->imgSize);
	}
	Matrix<float> KK(){
		return Matrix<float>(3, 3,  this->matrix );
	}
	CvMat *KK( CvMat *header ){
		return cvInitMatHeader( header, 3, 3, CV_32F, this->matrix );
	}
	Matrix<float> kc(){
		return Matrix<float>(4,1, this->distortion );
	}
	CvMat *kc( CvMat *header ){
		return cvInitMatHeader( header, 4, 1, CV_32F, this->distortion );
	}
	Matrix<float> R(){
		return Matrix<float>(3,3, this->rotMatr);
	}
	CvMat *R( CvMat *header ){
		return cvInitMatHeader( header, 3, 3, CV_32F, this->rotMatr );
	}
	Matrix<float> t(){
		return Matrix<float>(3,1,this->transVect);
	}
	CvMat *t( CvMat *header ){
		return cvInitMatHeader( header, 3, 1, CV_32F, this->transVect );
	}
};

struct StereoCamera: CvStereoCamera {
	StereoCamera(){
		camera[0] = (CvCamera *) cvAlloc( sizeof(CvCamera) );
		camera[1] = (CvCamera *) cvAlloc( sizeof(CvCamera) );
	}
	~StereoCamera(){
		cvFree(&camera[0]);
		cvFree(&camera[1]);
	}
	Camera & left(){
		return *((cv::Camera *) this->camera[0]);
	}
	Camera & right(){
		return *((cv::Camera *) this->camera[1]);
	}
	Matrix<float> F(){
		return Matrix<float>(3,3,  this->fundMatr);
	}
	CvMat *F( CvMat *header ){
		return cvInitMatHeader( header, 3, 3, CV_32F, this->fundMatr );
	}
	Matrix<float> R(){
		return Matrix<float>(3,3,  this->rotMatrix);
	}
	CvMat *R( CvMat *header ){
		return cvInitMatHeader( header, 3, 3, CV_32F, this->rotMatrix );
	}
	Matrix<float> t(){
		return Matrix<float>(3,1,  this->transVector);
	}
	CvMat *t( CvMat *header ){
		return cvInitMatHeader( header, 3, 1, CV_32F, this->transVector );
	}
	bool save( const char * filename ){
		int i, j;

		FILE* f = fopen( filename, "w" );

		if( !f ) return false;

		fprintf( f, "%d\n\n", 2 );

		for( i = 0; i < 2; i++ )
		{
			CvCamera * c = this->camera[i];
			for( j = 0; j < (int)(sizeof(*c)/sizeof(float)); j++ )
			{
				fprintf( f, "%15.10f ", ((float*)c)[j] );
			}
			fprintf( f, "\n\n" );
		}

		/* Save stereo params */

		/* Save quad */
		for( i = 0; i < 2; i++ )
		{
			for( j = 0; j < 4; j++ )
			{
				fprintf(f, "%15.10f ", this->quad[i][j].x );
				fprintf(f, "%15.10f ", this->quad[i][j].y );
			}
			fprintf(f, "\n");
		}

		/* Save coeffs */
		for( i = 0; i < 2; i++ )
		{
			for( j = 0; j < 9; j++ )
			{
				fprintf(f, "%15.10lf ", this->coeffs[i][j/3][j%3] );
			}
			fprintf(f, "\n");
		}


		fclose(f);
		return true;
	}
};

// output stream operators
std::ostream & operator<< (std::ostream & out, Camera & cam){
	out<<"imSize = ["<<cam.size().width<<","<<cam.size().height<<"]"<<std::endl;
	out<<"KK = "<<cam.KK()<<std::endl;
	out<<"kc = "<<cam.kc()<<std::endl;
	return out;
}
std::ostream & operator<< (std::ostream & out, StereoCamera & cam){
	out<<"Left Camera:\n"<<cam.left()<<std::endl;
	out<<"Right Camera:\n"<<cam.right()<<std::endl;
	out<<"R = "<<cam.R()<<std::endl;
	out<<"t = "<<cam.t()<<std::endl;
	return out;
}

} // namespace cv
