#ifndef __EOS_CC
#define __EOS_CC
#include "eos.h"

void EOS::freegascalc_onespecies_finitewidth(CboseMap &npidens,CboseMap &npiP,CboseMap &npiepsilon,CboseMap &npidedt,CparameterMap *parmap,CresInfo *resinfo,double T,double resmass,double m1, double m2, double reswidth_alpha,
                                              double &epsilon,double &P,double &dens,double &sigma2,double &dedt,double &maxweight){

    double E,E_S0,rho,rho_0,res_dens,weight,avg_weight,normal,N;
    double sum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
    int Ncounter = 0;
    double width=resinfo->width;
    double minmass=resinfo->minmass;
    maxweight=0.0;
    res_dens=gsl_sf_bessel_Kn(2,resmass/T)*resmass*resmass*T/(2*PI*PI*pow(HBARC,3.0));

    for (auto it = resinfo->spectmap.begin(); it != resinfo->spectmap.end(); it++) {
        E = (*it).first;
        rho = (*it).second;
        E_S0 = E - resmass;
        rho_0 = (1/PI)*(width/2.0)/(0.25*width*width+E_S0*E_S0);
        freegascalc_onespecies(npidens,npiP,npiepsilon,npidedt,parmap,resinfo,T,E,epsilon,P,dens,sigma2,dedt);
        weight=rho*dens/(rho_0*res_dens);
        if(weight>maxweight)  maxweight=weight;
        esum+=epsilon*rho/rho_0;
        psum+=P*rho/rho_0;
        dsum+=dens*rho/rho_0;
        sigsum+=sigma2*rho/rho_0;
        dedtsum+=dedt*rho/rho_0;
        sum+=rho/rho_0;
        Ncounter++;
    }

    N=Ncounter;
    avg_weight=sum/N;
    normal=1.0/avg_weight;
    epsilon=normal*esum/N;
    P=normal*psum/N;
    dens=normal*dsum/N;
    sigma2=normal*sigsum/N;
    dedt=normal*dedtsum/N;
}

void EOS::freegascalc_onespecies(CboseMap &npidens,CboseMap &npiP,CboseMap &npiepsilon,CboseMap &npidedt,CparameterMap *parmap,CresInfo *resinfo,double T,double m,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
    const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
    double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
    bool pion;
    int n;
    m2=m*m;
    m3=m2*m;
    m4=m2*m2;
    //this whole if/else loop does all the bose correction checks
    if (resinfo->code==-211 || resinfo->code==211 || resinfo->code==111) {
        pion=true;
        if (parmap->getB("BOSE_CORR",false)) {
            n=parmap->getI("N_BOSE_CORR",1);
        }
        else n=1;
    }
	else {
        n=1;
        pion=false;
    }

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
        P=epsilon=dens=dedt=0.0;
        for (int i=1;i<=n;i++) {
            double Ti=T/i;
            t2=Ti*Ti;
            t3=t2*Ti;
            z=m/Ti;
            k0=Bessel::K0(z);
            k1=Bessel::K1(z);
            //k0=boost::math::cyl_bessel_k(0.0,z);
            //k1=boost::math::cyl_bessel_k(1.0,z);

            double temp=prefactor*(m2*t2*k0+2.0*m*t3*k1);
            P+=temp;
            if(pion) npiP[resinfo->code].push_back(temp);

            dens+=temp/Ti;
            if(pion) npidens[resinfo->code].push_back(temp/Ti);

            temp=prefactor*(3.0*m2*t2*k0+(m3*Ti+6.0*m*t3)*k1);
            epsilon+=temp;
            if(pion) npiepsilon[resinfo->code].push_back(temp);

            k0prime=-k1;
            k1prime=-k0-k1/z;

            temp=prefactor*(6.0*m2*Ti*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/Ti)+6.0*m2*Ti)*k1prime);
            dedt+=temp;
            if(pion) npidedt[resinfo->code].push_back(temp);
        }
        z=m/T;
        Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
        I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
        I2=24.0*pow(T,5);
        sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
    }
}


#endif
