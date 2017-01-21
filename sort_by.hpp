#include <vector>
#include <algorithm>

/// sort_by_node  -- holds a copy of object and the value to be sorted by
template <typename obj_type, typename value_type>
struct sort_by_node {
	typedef sort_by_node<obj_type, value_type> self_type;

	obj_type obj;
	value_type val;
	
	sort_by_node(const obj_type &o, const value_type &v ):
		obj(o), val(v)
	{
	}
	
	inline bool operator< (const self_type & a) const{
		//std::cout<<*obj<<"="<<*val
		return this->val < a.val;
	}
	self_type & operator=(const self_type & n) {
		val = n.val;
		obj = n.obj;
		return (*this);
	}
};

/// stupid function to sort an array based on the values of a second array.
// could be more efficient about copying objects.  
template <typename iterator_T1, typename iterator_T2>
void sort_by( iterator_T1 first_obj, iterator_T1 last_obj,
		     iterator_T2 first_val, iterator_T2 last_val ){

	// copy elements into a new array which stores the object,value pair
	std::vector< sort_by_node<typename iterator_T1::value_type, typename iterator_T2::value_type> > array;

	int len_obj = last_obj - first_obj;
	int len_val = last_val - first_val;
	assert( len_val == len_obj );

	iterator_T1 obj = first_obj;
	iterator_T2 val = first_val;
	for(int i = 0; i<len_obj; i++){
		array.push_back( 
				sort_by_node<typename iterator_T1::value_type, typename iterator_T2::value_type >( *obj, *val ) ); 
		++val;
		++obj;
	}

	// now sort the array -- the '<' operator is defined in the sort_by_node struct 
	// to only compare the 'val' field
	std::sort( array.begin(), array.end() );

	// now copy the objects and values from the sorted array back to the original
	obj = first_obj;
	val = first_val;
	for(int i=0; i<len_obj; i++){
		*obj = array[i].obj;
		*val = array[i].val;
		++val;
		++obj;
	}
}
