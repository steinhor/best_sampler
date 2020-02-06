#ifndef __INCLUDE_GSLBESS_CC__
#define __INCLUDE_GSLBESS_CC__
#include "pratt_sampler/sf.h"

using namespace std;
using namespace pratt_sampler;

double Bessel::K0(double x){
	if(x>60.0)
		return 0.0;
	else
		return gsl_sf_bessel_K0(x);
}

double Bessel::K1(double x){
	if(x>60.0)
		return 0.0;
	else
		return gsl_sf_bessel_K1(x);
}

#endif
