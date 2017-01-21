#ifndef __SCALE_POINT_HH
#define __SCALE_POINT_HH
#include <cv/wincompat.h>
#include <cv/image.h>
#include <iostream>

struct ScalePoint {
	typedef cv::Image<float> image_t;

	float x;
	float y;
	float value;
	int level;
	int subLevel;
	float scale;
	float imgScale;
	float angle;
	
	ScalePoint(const int &X, const int &Y, const int &l, const int &sl, const float &s, const float &is) :
		x(X),y(Y),value(0),level(l),subLevel(sl),scale(s),imgScale(is) {
		}
	ScalePoint(const float &X, const float &Y, const int &l, const int &sl, const float &s, const float &is) :
		x(X),y(Y),value(0),level(l),subLevel(sl), scale(s), imgScale(is) {
		}
	int   imX()       const { return (int) (x*exp2(level-1)); }
	int   imY()       const { return (int) (y*exp2(level-1)); }
	float flX()       const { return x*exp2(level-1); }
	float flY()       const { return y*exp2(level-1); }
	
	
	void recalc(){
		imgScale = exp2(level-1)*scale;
	}
	
	// calculate radius of feature in pixel units in original image
	float getRadius() const { return ((((3.0 * imgScale) * 3.0) / 2.0) + 0.5)*.5; }
		
	friend std::ostream & operator << (std::ostream & out, const ScalePoint & p){
		//return out<<"["<<p.x<<","<<p.y<<":"<<p.level<<":"<<p.subLevel<<" v"<<p.value<<"]";
		//return out<<"["<<p.imX()<<", "<<p.imY()<<", "<<p.scale<<", "<<p.angle<<"]";
		return out<<"["<<(int)p.x<<", "<<(int)p.y<<", "<<p.subLevel<<", "<<p.level<<"]";
	}
	void calcGradientHistogram(float x, float y, float radius, const image_t & Ix, const image_t & Iy, double * hist, int nbins){
	
	}
	void draw(IplImage * im, const CvScalar &color, const char * label=NULL, CvFont * font=NULL){
		CvPoint p1 = cvPoint(this->imX(), this->imY());
		CvPoint p2 = cvPoint((int) (this->flX()+this->getRadius()*cos(this->angle)),
				(int)(this->flY()+this->getRadius()*sin(this->angle)));
		cvCircle(im, p1, (int) this->getRadius(), color);
		cvLine(im, p1, p2, color);
		if(label){
			cvPutText(im, label, p1, font, color);
		}
	}

	// calculate dominant orientation
	void calcOrientation(const image_t & Ix, const image_t & Iy){
		int nbins = 36;
		double hist[36];
		double sigma = this->scale * 3.0;
		double factor = 2.0*sigma*sigma;
		int radius = (int)(3.0*sigma / 2.0 + 0.5);
		//std::cout<<"radius = "<<radius<<std::endl;

		int radiusSq = radius*radius;

		bzero(hist, sizeof(double)*36);

		int xmin = cvFloor(MAX(0, this->x-radius));
		int ymin = cvFloor(MAX(0, this->y-radius));
		int xmax = cvFloor(MIN(Ix.width, this->x+radius)); 
		int ymax = cvFloor(MIN(Ix.height, this->y+radius));

		//std::cout<<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<<std::endl;
		for(int j=ymin; j<ymax; j++){
			for(int i=xmin; i<xmax; i++){
				int relX = cvRound(i - this->x);
				int relY = cvRound(j - this->y);

				// only consider points in circle
				double d = relX*relX + relY*relY;
				if(d>radiusSq) continue;

				// gaussian weight
				double gWeight = exp(  -( d / factor)  );

				int bin;
				float magnitude;
				float angle;
				calcGradient(i,j, Ix, Iy, magnitude, angle);
				bin = binGradient(angle, nbins);
				hist[bin] += magnitude*gWeight;
			}
		}
		averageWeakBins(hist, nbins);

		// find maximum bin
		double maxv=hist[0];
		int maxb=0;
		//std::cout<<"bin 0  = "<<hist[0]<<std::endl;
		for(int b=1; b<nbins; b++){
			//std::cout<<"bin "<<b<<" = "<<hist[b]<<std::endl;
			if(hist[b]>maxv){
				maxv=hist[b];
				maxb=b;
			}
		}
		this->angle = (M_PI*2*maxb)/nbins;
	}
	
	// Average the content of the direction bins.
	void averageWeakBins (double * hist, int nbins)
	{
		// TODO: make some tests what number of passes is the best. (its clear
		// one is not enough, as we may have something like
		// ( 0.4, 0.4, 0.3, 0.4, 0.4 ))
		for (int sn = 0 ; sn < 4 ; ++sn) {
			double firstE = hist[0];
			double last = hist[nbins - 1];

			for (int sw = 0 ; sw < nbins ; ++sw) {
				double cur = hist[sw];
				double next = (sw == (nbins - 1)) ?
					firstE : hist[(sw + 1) % nbins];

				hist[sw] = (last + cur + next) / 3.0;
				last = cur;
			}
		}
	}
	// Fit a parabol to the three points (-1.0 ; left), (0.0 ; middle) and
	// (1.0 ; right).
	//
	// Formulas:
	// f(x) = a (x - c)^2 + b
	//
	// c is the peak offset (where f'(x) is zero), b is the peak value.
	//
	// In case there is an error false is returned, otherwise a correction
	// value between [-1 ; 1] is returned in 'degreeCorrection', where -1
	// means the peak is located completely at the left vector, and -0.5 just
	// in the middle between left and middle and > 0 to the right side. In
	// 'peakValue' the maximum estimated peak value is stored.
	bool interpolateOrientation (double left, double middle,
			double right, double &degreeCorrection, double &peakValue)
	{
		double a = ((left + right) - 2.0 * middle) / 2.0;
		degreeCorrection = peakValue = 0;

		// Not a parabol
		if (a == 0.0)
			return (false);

		double c = (((left - middle) / a) - 1.0) / 2.0;
		double b = middle - c * c * a;

		if (c < -0.5 || c > 0.5) {
			std::cerr<<"ERROR: InterpolateOrientation: off peak ]-0.5 ; 0.5["<<std::endl;
			return false;
		}

		degreeCorrection = c;
		peakValue = b;

		return (true);
	}

	// TODO: FIXME!
	// First determine the real interpolated peak high at the maximum bin
	// position, which is guaranteed to be an absolute peak.
	//
	// XXX: should we use the estimated peak value as reference for the
	//   0.8 check or the original bin-value?
	void findAllOrientations(double * hist, bool * binIsKeypoint, int nbins, int maxBin, double peakRelThresh=0.8){
		double maxPeakValue, maxDegreeCorrection;
		interpolateOrientation (hist[maxBin == 0 ? (nbins - 1) : (maxBin - 1)],
				hist[maxBin], hist[(maxBin + 1) % nbins],
				maxDegreeCorrection, maxPeakValue);

		// Now that we know the maximum peak value, we can find other keypoint
		// orientations, which have to fulfill two criterias:
		//
		//  1. They must be a local peak themselves. Else we might add a very
		//     similar keypoint orientation twice (imagine for example the
		//     values: 0.4 1.0 0.8, if 1.0 is maximum peak, 0.8 is still added
		//     with the default threshhold, but the maximum peak orientation
		//     was already added).
		//  2. They must have at least peakRelThresh times the maximum peak
		//     value.
		for (int b = 0 ; b < nbins ; ++b) {
			binIsKeypoint[b] = false;

			// The maximum peak of course is
			if (b == maxBin) {
				binIsKeypoint[b] = true;
				continue;
			}

			// Local peaks are, too, in case they fulfill the threshhold
			if (hist[b] < (peakRelThresh * maxPeakValue))
				continue;

			int leftI = (b == 0) ? (nbins - 1) : (b - 1);
			int rightI = (b + 1) % nbins;
			if (hist[b] <= hist[leftI] || hist[b] <= hist[rightI])
				continue;   // no local peak

			binIsKeypoint[b] = true;
		}
	}

	template <typename InsertIterator>
	int calcOrientations(image_t & Ix, image_t & Iy, int width, int height, InsertIterator ii){
		int nbins = 36;
		double hist[36];
		bool binIsKeypoint[36];
		double sigma = this->scale * 3.0;
		double factor = 2.0*sigma*sigma;
		int radius = MIN((int)(3.0*sigma / 2.0 + 0.5), Ix.width/2);
		//std::cout<<"radius = "<<radius<<std::endl;

		int radiusSq = radius*radius;

		bzero(hist, sizeof(double)*36);

		// TODO -- these bounds seem wrong
		int xmin = (int) (MAX(1, this->x-radius) - this->x + Ix.width/2);
		int ymin = (int) (MAX(1, this->y-radius) - this->y + Ix.height/2);
		int xmax = (int) (MIN(width, this->x+radius) - this->x + Ix.width/2);
		int ymax = (int) (MIN(height, this->y+radius) - this->y + Ix.height/2);

		//std::cout<<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<<std::endl;
		for(int j=ymin; j<ymax; j++){
			for(int i=xmin; i<xmax; i++){
				int relX = i - Ix.width/2;
				int relY = j - Ix.height/2;
				//int imX = this->x + relX;
				//int imY = this->y + relY;
				// only consider points in circle
				double d = relX*relX + relY*relY;
				if(d>radiusSq) continue;

				// gaussian weight
				double gWeight = exp(  -( d / factor)  );

				int bin;
				float magnitude;
				float angle;
				calcGradient(i,j, Ix, Iy, magnitude, angle);
				bin = binGradient(angle, nbins);
				hist[bin] += magnitude*gWeight;
			}
		}

		averageWeakBins(hist, nbins);

		// find maximum bin
		double maxv=hist[0];
		int maxBin=0;
		//std::cout<<"bin 0  = "<<hist[0]<<std::endl;
		for(int b=1; b<nbins; b++){
			//std::cout<<"bin "<<b<<" = "<<hist[b]<<std::endl;
			if(hist[b]>maxv){
				maxv=hist[b];
				maxBin=b;
			}
		}
		
		// find all orientations
		findAllOrientations(hist, binIsKeypoint, nbins, maxBin);
		
		// now interpolate actual peak
		double oneBinRad = (2.0 * M_PI) / nbins;
		int numOrientations = 0;
		for(int b=0; b<nbins; b++){
			if(!binIsKeypoint[b]) continue;
			int bLeft = (b == 0) ? (nbins - 1) : (b - 1);
			int bRight = (b + 1) % nbins;

			// Get an interpolated peak direction and value guess.
			double peakValue;
			double degreeCorrection;

			if (interpolateOrientation (hist[bLeft], hist[b], hist[bRight],
						degreeCorrection, peakValue) == false)
			{
				//throw (new InvalidOperationException ("BUG: Parabola fitting broken"));
				std::cerr<<"BUG: Parabola fitting broken"<<std::endl;
				continue;
			}

			// [-1.0 ; 1.0] -> [0 ; binrange], and add the fixed absolute bin
			// position.
			// We subtract PI because bin 0 refers to 0, nbins-1 bin refers
			// to a bin just below 2PI, so -> [-PI ; PI]. Note that at this
			// point we determine the canonical descriptor anchor angle. It
			// does not matter where we set it relative to the peak degree,
			// but it has to be constant. Also, if the output of this
			// implementation is to be matched with other implementations it
			// must be the same constant angle (here: -PI).
			double degree = (b + degreeCorrection) * oneBinRad - M_PI;

			if (degree < -M_PI)
				degree += 2.0 * M_PI;
			else if (degree > M_PI)
				degree -= 2.0 * M_PI;

			ScalePoint s=(*this);
			s.angle = degree;
			*ii++ = s;
			numOrientations++;
		}
		return numOrientations;
	}

	static int binGradient(float angle, const int &nbins){
		while(angle<0) angle+=(2*M_PI);
		return (int)(nbins*angle*.5/M_PI);
	}
	static void calcGradient(const int &x, const int &y, const image_t & Ix, 
			const image_t & Iy, float & magnitude, float & angle){
		// compute gradient magnitude, orientation
		float dx = Ix[y][x];
		float dy = Iy[y][x];
		magnitude = sqrt( dx*dx + dy*dy );
		angle = atan2(dy,dx);
	}
	};

#endif //__SCALE_POINT_HH

