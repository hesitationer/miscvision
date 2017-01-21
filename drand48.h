#ifndef __DRAND48_H
#define __DRAND48_H

#include <stdlib.h>
#include <math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double drand48(){
	return ((double)rand())/(RAND_MAX+1.0);
}

// "Minimal" random number generator of Park and Miller with
// Bays-Durham shuffle and added safeguards. Returns a uniform random
// deviate between 0.0 and 1.0. Set or reset idum to any integer value
// (except the unlikely value MASK) to initialize the sequence; idum
// must not be altered between calls for successive deviates in a sequence. 
/*double ran0(int idum)
{
  long k;
  double ans;

  idum ^= MASK;           // XORing with MASK allows use of 0 and
  k = idum/IQ;            //     other simple bit patterns for idum.
  idum = IA * (idum-k*IQ) - IR*k; // Compute idum = (IA*idum) % IM without
  if (idum < 0) idum += IM;   //     overflows by Schrage's method.
  ans = AM * idum;        // Convert idum to a floating result.
  idum ^= MASK;           // Unmask before return.
  return ans;
}

//generate a gaussian sample - 
//courtesy of http://www-math.mit.edu/~spielman/ECC/randGen.c
static double gasdev()
{

        int iset=0;
        double gset;
        double fac,rsq,v1,v2;

        if (iset == 0) {
                do {
                        v1=2.0*drand48();
                        v2=2.0*drand48();
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq==0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        }
        else {
                iset =0;
                return gset;
        }
        return 0;
}*/
#endif // __DRAND48_H
