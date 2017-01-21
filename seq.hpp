#ifndef __CV_SEQ_H
#define __CV_SEQ_H

#include <cstddef>                // for ptrdiff_t, used next
#include <iterator>
#include <cxcore.h>

namespace cv{

template <typename T> struct SeqType{ static int type(){ return 0; } };
template <> struct SeqType<CvPoint>{ static int type(){ return CV_SEQ_ELTYPE_POINT; } };
template <> struct SeqType<uchar>{ static int type(){ return CV_SEQ_ELTYPE_CODE; } };
template <> struct SeqType<void *>{ static int type(){ return CV_SEQ_ELTYPE_PTR; } };
template <> struct SeqType<int>{ static int type(){ return CV_SEQ_ELTYPE_INDEX; } };
template <> struct SeqType<CvPoint2D32f>{ static int type(){ return CV_SEQ_ELTYPE_POINT; } };
template <> struct SeqType<CvPoint3D32f>{ static int type(){ return CV_SEQ_ELTYPE_POINT3D; } };


template <typename T, typename pointerT=T*> 
struct SeqIterator : public CvSeqReader{
	typedef int						        difference_type;
	typedef std::bidirectional_iterator_tag iterator_category;
	typedef T 								value_type;
	typedef pointerT						pointer;
	typedef T& 								reference;

	SeqIterator( ) : m_pos(0){
	}
	SeqIterator( CvSeq * seq ) : m_pos(0) {
		if(sizeof(T)!=seq->elem_size){
			OPENCV_ERROR( CV_StsUnmatchedFormats, __FUNCTION__, 
					"sizeof template argument is not equivalent to seq element size");
		}
		cvStartReadSeq( seq, this, 0);
	}
	SeqIterator( CvSeq * seq, int position ) : m_pos( position ){
		if(sizeof(T)!=seq->elem_size){
			OPENCV_ERROR( CV_StsUnmatchedFormats, __FUNCTION__, 
					"sizeof template argument is not equivalent to seq element size");
		}
		cvStartReadSeq( seq, this, 0);
		if(seq->total > 0) cvSetSeqReaderPos( this, position, 0);
	}
	reference operator*() const {
		return * (T*) this->ptr;
	}
	pointer operator->() const {
		return (T*) this->ptr;
	}
	SeqIterator<T> & operator++(){
		++m_pos;
		CV_NEXT_SEQ_ELEM( sizeof(T), (*this) );
		return (*this);
	}
	SeqIterator<T>  operator++( int ){
		SeqIterator<T> it = (*this);
		++(*this);
		return it;
	}
	SeqIterator<T>  operator--( int ){
		SeqIterator<T> it = (*this);
		--(*this);
		return it;
	}
	SeqIterator<T> & operator--(){
		--m_pos;
		CV_PREV_SEQ_ELEM( sizeof(T), (*this) );
		return (*this);
	}
	SeqIterator<T> & operator += (const difference_type & n){
		m_pos += n;
		cvSetSeqReaderPos( this, m_pos, 0);
		return (*this);
	}
	SeqIterator<T> & operator -= (const difference_type & n){
		m_pos -= n;
		cvSetSeqReaderPos( this, m_pos, 0);
	}
	SeqIterator<T> operator+( const difference_type & n ){
		return SeqIterator<T> ( *(this->seq), m_pos + n);
	}
	SeqIterator<T> operator-( const difference_type & n ){
		return SeqIterator<T> ( *(this->seq), m_pos - n );
	}
	difference_type operator-( const SeqIterator<T> & it ){
		return m_pos - it.m_pos;
	}
	bool operator!= ( const SeqIterator<T> & it ){
		return m_pos != it.m_pos;
	}
	bool operator == ( const SeqIterator<T> & it ){
		return m_pos == it.m_pos;
	}
	bool operator< ( const SeqIterator<T> & it ){
		return m_pos < it.m_pos;
	}
	bool operator> ( const SeqIterator<T> & it ){
		return m_pos > it.m_pos;
	}
	int m_pos;
};

template <typename T, typename pointerT=T*> 
class RevSeqIterator : public SeqIterator<T, pointerT> {
	RevSeqIterator( CvSeq * seq ){
		if(sizeof(T)!=seq->elem_size){
			OPENCV_ERROR( CV_StsUnmatchedFormats, __FUNCTION__, 
				"sizeof template argument is not equivalent to seq element size");
		}
		cvStartReadSeq( seq, this, 1 );
	}
};

template <typename T> 
class Seq : public CvSeq {
public:
	typedef SeqIterator<T, const T *> const_iterator;
	typedef RevSeqIterator<T, const T *> const_reverse_iterator;
	typedef SeqIterator<T> iterator;
	typedef RevSeqIterator<T> reverse_iterator;

	Seq(){
	}

	Seq(CvMemStorage * storage){
		this->init(storage);
	}
	void init(CvMemStorage * storage){
		// this is a little ugly...but no cvInitSeqHeader
		CvSeq * seq = cvCreateSeq(SeqType<T>::type(), sizeof(CvSeq), sizeof(T), storage);
		CvSeq * me = this;
		memcpy(me, seq, sizeof(CvSeq));
	}

	bool empty(){
		return this->total == 0;
	}
	~Seq(){
		this->clear();
	}
	T & front() {
		return (*this)[0];
	}
	const T & front() const{
		return (*this)[0];
	}
	T & back() {
		return (*this)[-1];
	}
	const T & back() const{
		return (*this)[-1];
	}
	void push_back(const T & elem){
		cvSeqPush(this, (void *) &elem);
	}
	void push(const T & elem){
		this->push_back( elem );
	}
	void pop_back(){
		T element;
		cvSeqPop(this, &element);
		return element;
	}
	void pop(const T & elem){
		this->pop_back( elem );
	}
	const T & operator [](int i) const {
		return *(T*)cvGetSeqElem(this, i);
	}
	T & operator [](int i){
		return *(T*)cvGetSeqElem(this, i);
	}
	void clear(){
		cvClearSeq(this);
	}
	int size() const{
		return this->total;
	}
	iterator begin() {
		return iterator( this );
	}
	iterator end() {
		return iterator( this, this->total );
	}
	const_iterator begin() const {
		return iterator( this );
	}
	const_iterator end() const {
		return iterator( this, this->total );
	}
};

} // namespace cv

#endif
