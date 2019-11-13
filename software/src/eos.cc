#ifndef __EOS_CC
#define __EOS_CC
#include "pratt_sampler/eos.h"

void EOS::freegascalc_onespecies_finitewidth(CresInfo *resinfo,double T,
double &epsilon,double &P,double &dens,double &sigma2,double &dedt,double &maxweight){
	double E,rho,weight,k2,k2mr;
	double rhosum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
	int n;
	maxweight=0.0;
	k2mr=resinfo->mass*resinfo->mass*Bessel::Kn(2,resinfo->mass/T);
	for (n=0;n<CresInfo::NSPECTRAL;n++){
		E=resinfo->GetEofN(n);
		if(E>0.0){
			rho=resinfo->spectvec[n];
			freegascalc_onespecies(T,E,epsilon,P,dens,sigma2,dedt);
			k2=E*E*Bessel::Kn(2,E/T);
			weight=rho*k2/k2mr;
			if(weight>maxweight){
				maxweight=1.0001*weight;  // safety factor of 1.0001
			}
			esum+=epsilon*rho;
			psum+=P*rho;
			dsum+=dens*rho;
			sigsum+=sigma2*rho;
			dedtsum+=dedt*rho;
			rhosum+=rho;
		}
	}
	epsilon=esum/rhosum;
	P=psum/rhosum;
	dens=dsum/rhosum;
	sigma2=sigsum/rhosum;
	dedt=dedtsum/rhosum;
}

void EOS::freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3;//I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	z=m/T;
	if(z>1000.0){
		P=epsilon=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,T);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
			exit(1);
		}
		t2=T*T;
		t3=t2*T;
		z=m/T;
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		dens=P/T;
		epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
	}
	/*
	z=m/T;
	Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
	I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
	I2=24.0*pow(T,5);
	sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	//printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
	*/
}


#endif
