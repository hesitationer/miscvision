#include <iostream>
#include <vector>
#include <kdtree++/kdtree.hpp>

class point_struct {
public:
	typedef subtype value_type;
	const inline value_type & operator[](const int N) const { return (*d)[N]; }
	t *d;
	int idx;
};

template<typename data_t, typename subtype>
inline bool operator==(point_struct<data_t, subtype> const& A, point_struct<data_t, subtype> const& B) {
  return A[0] == B[0];
}

template<typename data_t, typename subtype>
std::ostream& operator<<(std::ostream& out, point_struct<data_t, subtype> const& T)
{
    return out << '(' << T[0] << ')';
}

/*template <typename descriptor_t>
static void match_scale_points(iterator & A, iterator & B, std::vector<ScalePoint *> match, const double & epsilon){
	typedef point_struct<point_t, subtype> ps_t;

	KDTree::KDTree<__size, ps_t > tree;
	// insert everything in A to the kd-tree
	int size = A.size();
	for(int i=0; i<size; i++){
		ps_t p;
		p.d = &(A[i]);
		p.idx=i;
		tree.insert(p);
	}

	//tree.optimize();

	size = B.size();
	typename KDTree::KDTree<1, ps_t >::const_iterator it;
	for(int i=0; i<size; i++){
		// find the closest point
		ps_t p;
		std::vector<ps_t> v;
		p.d = &(B[i]);
		tree.find_within_range(p, 1, std::back_inserter(v));
		typename std::vector< ps_t >::const_iterator ci = v.begin();
		if (ci != v.end()){
			idx[i] = ci->idx;
			tree.erase(*ci);
		}
		else{
			idx[i] = -1;
		}
	}

}*/

// input two lists of features, a maximum distance to search
// output -- vector that, for each feature in list A, what index does it correspond to in list B, and its distance 

template <typename point_t, typename subtype, int __size>
class match_points {
public:
	static void match(std::vector<point_t> & A, std::vector<point_t> & B, std::vector<int> & idx, const subtype & epsilon){
		typedef point_struct<point_t, subtype> ps_t;

		KDTree::KDTree<__size, ps_t > tree;
		// insert everything in A to the kd-tree
		int size = A.size();
		for(int i=0; i<size; i++){
			ps_t p;
			p.d = &(A[i]);
			p.idx=i;
			tree.insert(p);
		}

		//tree.optimize();

		size = B.size();
		typename KDTree::KDTree<1, ps_t >::const_iterator it; 
		for(int i=0; i<size; i++){
			// find the closest point
			ps_t p;
			std::vector<ps_t> v;
			p.d = &(B[i]);
			tree.find_within_range(p, 1, std::back_inserter(v));
			typename std::vector< ps_t >::const_iterator ci = v.begin();
			if (ci != v.end()){
				idx[i] = ci->idx;
				tree.erase(*ci);
			}
			else{
				idx[i] = -1;
			}
		}
	}
};
