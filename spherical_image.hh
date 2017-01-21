template <typename T> class Hemisphere :: Image<T> {
	int _length;
	/* constructor
	 * @param n is the length in pixels of the diameter of the sphere
	 */
	Hemisphere(int n) : Image<T>(n/4,n), _length(n)
	{}

	/** accessor for point on surface of sphere
	 * @param theta azimuthal angle in x-y plane from x-axis (0->2PI)
	 * @param phi polar angle from the z-axis (0->PI)
	 */
	T & operator (float R, float theta) {
			
	}
	/** i * PI/N = phi 
	    rlength =~ # of pixels for this row
		rlength =~ diameter of sphere at phi
		        = sin (phi) * N
	*/
	int rlength(int i){
		return (int) (sin((M_PI * i)/ _length) *_length) + 1;
	}
}
