/** 
 * @brief Return the nearest neighbor to p in the range __first, __last
 * @param p Reference point to find nearest neighbor to
 * @param __first An input iterator indicating start of range 
 * @param __last An input iterator indicating the end of range
 */
template <typename point_t, typename _InputIter>
_InputIter find_nearest(const point_t & p, _InputIter __first, _InputIter __last ){
	double mindist=dist(p,*__first);
	_InputIter minit = __first;
	while(++__first!=__last){
		double d = dist(p, *__first);
		if(d<mindist){
			minit = __first;
			mindist = d;
		}
	}
	return minit;
}

template <typename point_t>
int find_nearest(const point_t &p, const std::vector<point_t> & v, double epsilon, double * rmin){
	double mindist=epsilon;
	int idx = -1;
	for(int i=0; (size_t)i<v.size(); i++){
		double d = dist(p, v[i]);
		//std::cerr<<"dist:"<<d<<std::endl;
		//std::cout<<p<<std::endl;
		//std::cout<<v[i]<<std::endl;
		if(d<mindist){
			idx=i;
			mindist=d;
		}
	}
	if(rmin) *rmin = mindist;
	return idx;
}
