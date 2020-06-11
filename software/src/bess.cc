#ifndef __INCLUDE_GSLBESS_CC__
#define __INCLUDE_GSLBESS_CC__
#include "msu_sampler/sf.h"

using namespace msu_sampler;

double Bessel::K0(double x){
	if(x<0 || x!=x){
		printf("K0(x<0), x=%g\n",x);
		exit(1);
	}
	if(x>50.0)
		return 0.0;
	else
		return gsl_sf_bessel_K0(x);
}

double Bessel::K1(double x){
	if(x<0 || x!=x){
		printf("K1(x<0), x=%g\n",x);
		exit(1);
	}
	if(x>50.0)
		return 0.0;
	else
		return gsl_sf_bessel_K1(x);
}

double Bessel::Kn(int n,double x){
	if(x<0 || x!=x){
		printf("Kn(x<0), n=%d, x=%g\n",n,x);
		exit(1);
	}
	if(x>50.0)
		return 0.0;
	else
		return gsl_sf_bessel_Kn(n,x);
}

#endif
