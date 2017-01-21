#ifndef __FIXED_PRIORITY_QUEUE_HPP
#define __FIXED_PRIORITY_QUEUE_HPP

#include <queue>
#include <vector>

/// simple wrapper to make data type weighted
template <typename data_t, typename weight_t=float>
struct weighted_var {
	data_t data;
	weight_t weight;

	weighted_var( const data_t & v, const weight_t & w):
		data(v), weight(w)
	{
	}
	weighted_var( const data_t & v ):
		data(v)
	{
	}
	data_t & operator ()() {
		return this->data;
	}
};

/// comparison function for weighted data type
template <typename data_t, typename weight_t>
struct compare_weighted_var {
	typedef weighted_var<data_t, weight_t> first_argument_type;
	typedef weighted_var<data_t, weight_t> second_argument_type;
	bool operator () ( const weighted_var<data_t, weight_t> & x1, const weighted_var<data_t, weight_t> & x2 ) const {
		return x1.weight < x2.weight;
	}
};

template <typename binary_logical>
struct binary_logical_not{
	typedef typename binary_logical::first_argument_type first_argument_type;
	typedef typename binary_logical::second_argument_type second_argument_type;
	binary_logical func;
	bool operator () ( const first_argument_type & x1, const second_argument_type & x2 ){
		return ! func( x1, x2 );
	}
};

/// simple priority queue modified to hold a fixed number of items
// note that items in pq are reversed from expected -- minimum item is the top
// and maximum is a leaf
template <typename T, class Sequence=std::vector< weighted_var<T, float> >, class Compare=compare_weighted_var<T, float> >
class fixed_priority_queue : public std::priority_queue <T, Sequence, binary_logical_not<Compare> > {
	size_t maxSize;
public:
	typedef std::priority_queue <T, Sequence, binary_logical_not<Compare> > base_type;
	typedef T value_type;
	typedef weighted_var<T,float> weighted_type;
	typedef size_t size_type;
	typedef typename Sequence::iterator iterator;
	typedef typename Sequence::const_iterator const_iterator;

	fixed_priority_queue( size_t sz ):
		base_type(),
		maxSize( sz )
	{
	}
	void push( const value_type & v, float weight ){
		this->push( weighted_type( v, weight ) );
	}
	void push( const weighted_type & v ){
		if( this->size() >= maxSize ){
			if(!this->comp( this->top(), v ) ){
				this->pop();
			}
			else{
				return;
			}
		}
		base_type::push( v );
	}

	typename Sequence::value_type & back(){
		return this->c.back();
	}
	const typename Sequence::value_type & back() const{
		return this->c.back();
	}
	iterator begin() {
		return this->c.begin();
	}
	iterator end() {
		return this->c.end();
	}
	const_iterator begin() const{
		return this->c.begin();
	}
	const_iterator end() const{
		return this->c.end();
	}
	void clear(){
		this->c.clear();
	}
	void sort(){
		std::sort_heap( this->begin(), this->end(), this->comp );
	}
	void reserve( size_t size ){
		maxSize = size;
	}
};

#endif // __FIXED_PRIORITY_QUEUE_HPP
