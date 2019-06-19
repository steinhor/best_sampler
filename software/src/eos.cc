#ifndef __EOS_CC
#define __EOS_CC
#include "eos.h"



void EOS::freegascalc_onespecies_finitewidth(double T,double resmass, double m1, double m2, double width, double reswidth_alpha,double spin_deg,
                                              double minmass,double &epsilon,double &P,double &dens,double &sigma2,double &dedt,double &maxweight){
    
    double kr,k,E0,gamma,rho,rho_0,gammas,n0,res_dens,weight,avg_weight,normal;
    double sum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
    int N = 100;
    // E_S0 is an array of E' values where E'=E-resmass
    double E_S0;
    double E;
    
    // E0 = minmass if minmass declared already in resonance.cc
    E0=m1+m2;
    maxweight=0.0;
    res_dens=gsl_sf_bessel_Kn(2,resmass/T)*resmass*resmass*T/(2*PI*PI*pow(HBARC,3.0));
    kr=sqrt(abs(pow((resmass*resmass-m1*m1-m2*m2),2.0)-4.0*m1*m1*m2*m2))/(2.0*resmass);
    //printf("kr=%g\n",kr);
    //printf("E0=%g\n",E0);
    for(int n=0;n<N;n++)
    {
        double Sum_E=(n+0.5)/N;
        E_S0 = 0.5*width*tan(PI*(Sum_E - .5));
        if((E_S0+resmass)>=minmass)
        {
            E = E_S0+resmass;
            //printf("n=%d\n",n);
            //printf("E=%g\n",E);
            
            
            
            k=sqrt(abs(pow((E*E-m1*m1-m2*m2),2.0)-(4.0*m1*m1*m2*m2)))/(2.0*E);
           
            if(spin_deg<1.001)
            {
                gamma=width*(resmass/E)*(k/kr);
            }
            else
            {
                gamma=width*(resmass/E)*((k*k*k)/(kr*kr*kr))*((kr*kr+HBARC*HBARC)/(k*k+HBARC*HBARC));
            }
            
            
            rho=(2.0)/(width*PI)*0.25*gamma*gamma/((0.25*gamma*gamma)+(resmass-E)*(resmass-E));
            
            
            
            rho_0 = (1/PI)*(width/2.0)/(0.25*width*width+E_S0*E_S0);
            
            freegascalc_onespecies(T,E,epsilon,P,dens,sigma2,dedt);
            
            //printf("weight=%g\n",weight);
            //printf("rho_0=%g\n",rho_0);
            
            //printf("T=%g,E=%g,epsilon=%g,P=%g,dens=%g,sigma2=%g,dedt=%g\n",T,E,epsilon,P,dens,sigma2,dedt);
            
            weight=rho*dens/(rho_0*res_dens);
           
            if(weight>maxweight)
                maxweight=weight;
            //printf("dens=%g\n",dens);
            //printf("rho=%g\n",rho);
            //printf("rho_0=%g\n",rho_0);
            esum+=epsilon*rho/rho_0;
            psum+=P*rho/rho_0;
            dsum+=dens*rho/rho_0;
            //printf("dsum=%g\n",dsum);
            sigsum+=sigma2*rho/rho_0;
            dedtsum+=dedt*rho/rho_0;
            sum+=rho/rho_0;
            //printf("sum=%g\n",sum);
        }
        else{
            // blank
        }
    }
    
    avg_weight=sum/N;
    normal=1.0/avg_weight;
    epsilon=normal*esum/N;
    P=normal*psum/N;
    dens=normal*dsum/N;
    //printf("dens=%g\n",dens);
    sigma2=normal*sigsum/N;
    dedt=normal*dedtsum/N;
    /*
    epsilon=esum/N;
    P=psum/N;
    dens=dsum/N;
    //printf("dens=%g\n",dens);
    sigma2=sigsum/N;
    dedt=dedtsum/N;
    */
}

void EOS::freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
    const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
    double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
    m2=m*m;
    m3=m2*m;
    m4=m2*m2;
    t2=T*T;
    t3=t2*T;
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
        k0=Bessel::K0(z);
        k1=Bessel::K1(z);
        //k0=boost::math::cyl_bessel_k(0.0,z);
        //k1=boost::math::cyl_bessel_k(1.0,z);
        P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
        epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
        dens=P/T;
        k0prime=-k1;
        k1prime=-k0-k1/z;
        dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
        Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
        I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
        I2=24.0*pow(T,5);
        sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
        //printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
    }
}


#endif
