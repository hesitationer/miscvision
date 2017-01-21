struct Mat5Header {
	char description[124];
	short version;
	short endianess;
};

// borrowed from octave
enum Mat5Type
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


union Mat5SmallElement {
	struct {
		short type;
		short size;
	} 
	int type;
	char bytes[4];
};
	
class   GrFmtMatReader : public GrFmtReader
{
public:

    GrFmtMatReader( const char* filename );
    virtual ~GrFmtMatReader();

    virtual bool  ReadHeader();
    virtual bool  ReadData( uchar* data, int step, int color );
    virtual void  Close();

protected:
	RLByteStream * m_stream;
	bool m_swap;
};

bool GrFmtMatReader::GrFmtMatReader( const char * filename ) : GrFmtReader( filename ) {
}

bool GrFmtMatReader::ReadHeader(){
	Mat5Header header;
	Mat5SmallElement element;

    assert( strlen(m_filename) != 0 );
	m_stream = new RLByteStream;
    if( !m_stream->Open( m_filename )) return false;

	// read in header
	m_stream->ReadBytes(header, sizeof(Mat5Header)); 
	
	// check endianness
	// may have to open an RMByteStream
	if (header.endianess == 0x4d49){
		m_swap = false;
	}
	else if (header.endianess == 0x494d){
		m_swap = true;
	}
	else
	{
		fprintf(stderr, "load: can't read binary file");
		return false;
    }

	// fix version number
	if (! m_swap){
		// version number is inverse m_swapped!
		version = ((version >> 8) & 0xff) + ((version & 0xff) << 8);
	}
	
	if (version != 1){
		warning ("load: found version %d binary MAT file, "
				"but only prepared for version 1", version);	
	}

	// now read element header
	this->ReadTag( &type, &element_length );

	if(type != miMATRIX){
		fprintf(stderr, "Invalid element type -- expected miMATRIX, found %d\n", type);
		return false;
	}
	
	// empty array
	if(element_length==0){
		m_width=0;
		m_height=0;
		return true;
	}

	// read array sub element 
	this->ReadTag( &type, &len );
	if
	return true;
}

bool ReadTag( int * type, int * len ){
	// first check if it is a small element
	int tmp;
	int upper;
	if(!m_stream->ReadBytes(&tmp, 4)){
		return false;
	}
	
	upper = (tmp >> 16) & 0xffff;
	*type = temp & 0xffff;
	if(upper){
		*len = upper;
	}
	else{
		if(!m_stream->ReadBytes(&tmp, 4)){
			return false;
		}
		*len = tmp;
	}
	return true;	
}

bool ReadData( uchar * data, int step, int color ){
	
}

void GrFmtMatReader::Close(){
	m_stream->Close();
}
