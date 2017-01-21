#include <stdio.h>
#include <cxcore.h>
#include <cxmisc.h>

typedef struct Mat5Header {
    char description[116];
	char data_offset[8];
    short version;
    char endianess[2];
} Mat5Header; 

// structs borrowed from octave
enum Mat5ArrayType
  {
    mxCELL_CLASS=1,     // cell array
    mxSTRUCT_CLASS,     // structure
    mxOBJECT_CLASS,     // object
    mxCHAR_CLASS,       // character array
    mxSPARSE_CLASS,     // sparse array
    mxDOUBLE_CLASS,     // double precision array
    mxSINGLE_CLASS,     // single precision floating point
    mxINT8_CLASS,       // 8 bit signed integer
    mxUINT8_CLASS,      // 8 bit unsigned integer
    mxINT16_CLASS,      // 16 bit signed integer
    mxUINT16_CLASS,     // 16 bit unsigned integer
    mxINT32_CLASS,      // 32 bit signed integer
    mxUINT32_CLASS,     // 32 bit unsigned integer
    mxINT64_CLASS,      // 64 bit signed integer
    mxUINT64_CLASS,     // 64 bit unsigned integer
    mxFUNCTION_CLASS            // Function handle
  };
enum Mat5DataType
  {
    miINT8 = 1,         // 8 bit signed
    miUINT8,            // 8 bit unsigned
    miINT16,            // 16 bit signed
    miUINT16,           // 16 bit unsigned
    miINT32,            // 32 bit signed
    miUINT32,           // 32 bit unsigned
    miSINGLE,           // IEEE 754 single precision float
    miRESERVE1,
    miDOUBLE,           // IEEE 754 double precision float
    miRESERVE2,
    miRESERVE3,
    miINT64,            // 64 bit signed
    miUINT64,           // 64 bit unsigned
    miMATRIX            // MATLAB array
  };


typedef struct Mat5Element {
    int data_type;
	int size;
} Mat5Element;

CV_IMPL int cvSaveMat5(const char * fname, const CvArr * arr){
	int ret = 0;
	CV_FUNCNAME( "cvSaveMat5" );
	
	__BEGIN__;

	CvMat * mat = (CvMat *) arr;
	CvMat matStub;
	Mat5Header header;
	Mat5Element element;
	Mat5Element flagselement;
	Mat5Element dimelement;
	Mat5Element nameelement;
	Mat5Element dataelement;
	FILE * f;
	int ndim = 2;
	int elemsize;
	int array_flags;

	f = fopen(fname, "wb");
	if( !f ){
		CV_ERROR( CV_StsError, "Couldn't open file for writing");
		return -1;
	}

	if( !CV_IS_MAT( mat ) ){
		mat = cvGetMat(mat, &matStub);
	}

	// configure header
	snprintf(header.description, sizeof(header.description), "MATLAB 5.0 MAT-file, Created by: OpenCV");
	memset(header.data_offset, 0, sizeof(header.data_offset));
	header.version = 0x0100;  // as per spec 
	header.endianess[0] = 'I';
	header.endianess[1] = 'M';

	// configure element header
	element.data_type = miMATRIX;
	element.size = 0;

	// configure array flags 
	flagselement.data_type = miUINT32;
	flagselement.size = 8;
	element.size += sizeof(flagselement) + flagselement.size;

	switch( CV_MAT_DEPTH(mat->type) ){
	case CV_8S:
		array_flags = mxINT8_CLASS;
		dataelement.data_type = miINT8;
		break;
	case CV_8U:
		array_flags = mxUINT8_CLASS;
		dataelement.data_type = miUINT8;
		break;
	case CV_16S:
		array_flags = mxINT16_CLASS;
		dataelement.data_type = miINT16;
		break;
	case CV_16U:
		array_flags = mxUINT16_CLASS;
		dataelement.data_type = miUINT16;
		break;
	case CV_32S:
		array_flags = mxINT32_CLASS;
		dataelement.data_type = miINT32;
		break;
	case CV_32F:
		array_flags = mxSINGLE_CLASS;
		dataelement.data_type = miSINGLE;
		break;
	case CV_64F:
		array_flags = mxDOUBLE_CLASS;
		dataelement.data_type = miDOUBLE;
		break;
	default:
		CV_ERROR( CV_StsUnsupportedFormat, "" );
	}
	
	// dimensions sub-element
	if(CV_MAT_CN( mat->type ) > 1){
		ndim = 3;
	}
	dimelement.data_type = miINT32;
	dimelement.size = ndim*4;
	element.size += sizeof(dimelement) + dimelement.size;

	// name sub-element
	nameelement.data_type = miINT8;
	nameelement.size = 8;

	element.size += sizeof(nameelement) + nameelement.size;

	// data sub-element
	dataelement.size = mat->rows*mat->cols*CV_ELEM_SIZE(mat->type);
	element.size += dataelement.size + sizeof(dataelement);
	
	// write header
	fwrite( &header, sizeof(header), 1, f);

	// write element header
	fwrite( &element, sizeof(element), 1, f);
	
	// array flags element header
	fwrite( &flagselement, sizeof(flagselement), 1, f);

	// array flags (8 bytes, second 4 bytes undefined, just write twice)
	fwrite( &array_flags, sizeof(array_flags), 1, f);
	fwrite( &array_flags, sizeof(array_flags), 1, f);

	// dimensions header
	fwrite( &dimelement, sizeof(dimelement), 1, f);
	
	// dimension array
	fwrite( &(mat->rows), sizeof(int), 1, f);
	fwrite( &(mat->cols), sizeof(int), 1, f);
	if(ndim>2){
		int cn = CV_MAT_CN( mat->type );
		fwrite( &cn, sizeof(int), 1, f);
	}

	// name header
	fwrite( &nameelement, sizeof(nameelement), 1, f);
	char arrayname[8];
	sprintf( arrayname, "CVMAT" );
	fwrite( arrayname, sizeof(char), 8, f);

	// data header
	fwrite( &dataelement, sizeof(dataelement), 1, f);
	
	// actual data -- stupid matlab packs arrays opposite of opencv
	elemsize=CV_ELEM_SIZE(mat->type);
	for(int j=0; j<mat->cols; j++){
		for(int i=0; i<mat->rows; i++){
			fwrite( mat->data.ptr+i*mat->step+elemsize*j, elemsize, 1, f );
		}
	}

	fclose(f);

	__END__;
	return ret;
}

#ifdef MAIN
#include <cv/matrix.hh>
int main(int argc, char ** argv){
	Matrix<float> mat(3, 14);
	for(int i=0; i<mat.rows; i++){
		for(int j=0; j<mat.cols; j++){
			mat[i][j] = drand48();
		}
	}
	std::cout<<mat<<std::endl;
	cvSaveMat5("A.mat", &mat);
}
#endif
