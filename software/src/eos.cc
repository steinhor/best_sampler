#ifndef __EOS_CC
#define __EOS_CC
#include "msu_sampler/eos.h"
#include "msu_sampler/sf.h"
#include "msu_sampler/constants.h"

using namespace msu_sampler;

void EOS::freegascalc_onespecies_finitewidth(CresInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &dedt,double &p4overE3){
	int iE,nE;
	double E,rho;
	double rhosum=0.0,esum=0.0,psum=0.0,dsum=0.0,p4overE3sum=0.0,dedtsum=0.0;
	nE=resinfo->SpectVec.size();
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		rho=resinfo->SpectVec[iE];
		freegascalc_onespecies(T,E,epsilon,P,dens,dedt,p4overE3);
		esum+=epsilon*rho;
		psum+=P*rho;
		dsum+=dens*rho;
		dedtsum+=dedt*rho;
		p4overE3sum+=p4overE3;
		rhosum+=rho;
	}
	epsilon=esum/rhosum;
	P=psum/rhosum;
	dens=dsum/rhosum;
	dedt=dedtsum/rhosum;
	p4overE3=p4overE3sum/rhosum;
}

void EOS::freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &dedt,double &p4overE3){
	const double prefactor=1.0/(2.0*M_PI*M_PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	z=m/T;
	if(z>40){ // non-relativistic approximations
		dens=exp(-z)*pow(m*T/(2.0*M_PI*HBARC*HBARC),1.5);
		P=dens*T;
		epsilon=(m+1.5*T)*dens;
		dedt=z*z*dens+3.0*z*dens;
		p4overE3=15.0*dens*T*T/m;
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
	p4overE3=Getp4overE3(T,m);
}

double EOS::Getp4overE3(double T,double m){
	double z,J=0.0,nfact,sign;
	double p4overE3=0.0;
	const int nmax=50;
	double G[nmax+5];
	int n;
	z=m/T;
	if(z<0 || z!=z){
		printf("K0(z<0)\n, z=%g\n",z);
		exit(1);
	}
	if(z<50.0){
		G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
		for(int j=1;j<nmax+5;j++){
			n=5-2*j;
			if(n!=-1)
				G[j]=(-exp(-z)/n)+(G[j-1]*z*z-z*exp(-z))/((n+1.0)*n);
			else G[j]=gsl_sf_gamma_inc(-1,z)*z;
		}
		J=0.0;
		nfact=1.0;
		sign=1.0;
		for(n=0;n<nmax;n+=1){
			if(n>0)
				sign=-1.0;
			J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
			nfact=nfact*0.5/(n+1.0);
			if(n>0)
				nfact*=(2.0*n-1.0);
		}
		//p4overE3=pow(m,4)*(-z*J+15.0*Bessel::Kn(2,z)/(z*z));
		p4overE3=pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
		p4overE3=2.0*p4overE3/(4.0*M_PI*M_PI*HBARC*HBARC*HBARC);
	}
	else
		p4overE3=0.0;
	
	/*
	double delp=0.002,pmag,e;
	double p4overE3test;
	printf("__ TESTING, m=%g, p4overE3=%g, J=%g, T=%g, z=%g\n",m,p4overE3,J,T,z);
	p4overE3test=0.0;
	for(pmag=0.5*delp;pmag<10.5;pmag+=delp){
		e=sqrt(m*m+pmag*pmag);
		//p4overE3test+=degen*(4.0*M_PI/(8.0*M_PI*M_PI*M_PI*HBARC*HBARC*HBARC))*pmag*pmag*dp*exp(-e/T)*( (2.0/3.0)*(pmag*pmag/e) - (2.0/15.0)*pow(pmag,4)/pow(e,3) );
		p4overE3test+=(4.0*M_PI/(8.0*M_PI*M_PI*M_PI*HBARC*HBARC*HBARC))*pmag*pmag*delp*exp(-e/T)*pow(pmag,4)/pow(e,3);
	}
	p4overE3=p4overE3test;
	//printf("XXXXXXX p4overE3test=%g, p4overE3test/p4overE3=%g, z=%g\n",p4overE3test,p4overE3test/p4overE3,z);
	*/
		
	return p4overE3;
}

#endif
