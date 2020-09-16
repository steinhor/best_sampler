#ifndef __INCLUDE_GSLBESS_CC__
#define __INCLUDE_GSLBESS_CC__
#include <cmath>
#include "msu_sampler/sf.h"
#include "msu_sampler/constants.h"

using namespace msu_sampler;
using namespace std;

double Bessel::K0(double x){
	if(x<0 || x!=x){
		printf("K0(x<0), x=%g\n",x);
		exit(1);
	}
	if(x>30.0){
		double answer=0.0,ak=1.0,factor;
		const int kmax=10;
		const double PI=4.0*atan(1.0);
		int k;
		for(k=0;k<=kmax;k++){
			if(k>0){
				factor=-(2.0*k-1)*(2.0*k-1);
				ak=ak*factor/(8.0*k);
			}
			answer+=ak/pow(x,k);
		}
		answer*=sqrt(0.5*PI/x)*exp(-x);
		//printf("K_0(%g)= %g =? %g, ratio=%g\n",x,gsl_sf_bessel_K0(x),answer,answer/gsl_sf_bessel_K0(x));
		return answer;
	}
	else
		return gsl_sf_bessel_K0(x);
}

double Bessel::K1(double x){
	if(x<0 || x!=x){
		printf("K1(x<0), x=%g\n",x);
		exit(1);
	}
	if(x>30.0){
		double answer=0.0,ak=1.0,factor;
		const int kmax=10;
		const double PI=4.0*atan(1.0);
		int k;
		for(k=0;k<=kmax;k++){
			if(k>0){
				factor=4.0-(2.0*k-1)*(2.0*k-1);
				ak=ak*factor/(8.0*k);
			}
			answer+=ak/pow(x,k);
		}
		answer*=sqrt(0.5*PI/x)*exp(-x);
		//printf("K_1(%g)= %g =? %g, ratio=%g\n",x,gsl_sf_bessel_K1(x),answer,answer/gsl_sf_bessel_K1(x));
		return answer;
	}
	else
		return gsl_sf_bessel_K1(x);
}

double Bessel::Kn(int n,double x){
	if(x<0 || x!=x){
		printf("Kn(x<0), n=%d, x=%g\n",n,x);
		exit(1);
	}
	if(x>30.0){
		double answer=0.0,ak=1.0,factor;
		const int kmax=10;
		const double PI=4.0*atan(1.0);
		int k;
		for(k=0;k<=kmax;k++){
			if(k>0){
				factor=4.0*n*n-(2.0*k-1)*(2.0*k-1);
				ak=ak*factor/(8.0*k);
			}
			answer+=ak/pow(x,k);
		}
		answer*=sqrt(0.5*PI/x)*exp(-x);
		//printf("K_%d(%g)= %g =? %g, ratio=%g\n",n,x,gsl_sf_bessel_Kn(n,x),answer,answer/gsl_sf_bessel_Kn(n,x));
		return answer;
	}
	else
		return gsl_sf_bessel_Kn(n,x);
}

#endif
