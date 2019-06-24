#ifndef __EOS_CC
#define __EOS_CC
#include "eos.h"




void EOS::freegascalc_onespecies_finitewidth(CresInfo *resinfo,double T,double resmass, double m1, double m2, double width, double reswidth_alpha,double spin_deg,
                                              double minmass,double &epsilon,double &P,double &dens,double &sigma2,double &dedt,double &maxweight){
    
    double kr,k,E0,rho,rho_0,gammas,n0,res_dens,weight,avg_weight,normal;
    double sum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
    double gamma=0;
    int N = 100;
    // E_S0 are E' values where E'=E-resmass
    double E_S0;
    double E;
    
    resmass=resinfo->mass;
    width=resinfo->width;
    minmass=resinfo->minmass;
    m1=resinfo->branchlist[0]->resinfo[0]->mass;
    m2=resinfo->branchlist[0]->resinfo[1]->mass;
    
    // E0 = minmass if minmass declared already in resonance.cc
    E0=m1+m2;
    maxweight=0.0;
    res_dens=gsl_sf_bessel_Kn(2,resmass/T)*resmass*resmass*T/(2*PI*PI*pow(HBARC,3.0));
    kr=sqrt(abs(pow((resmass*resmass-m1*m1-m2*m2),2.0)-4.0*m1*m1*m2*m2))/(2.0*resmass);
    
    for(int n=0;n<N;n++)
    {
        double Sum_E=(n+0.5)/N;
        E_S0 = 0.5*width*tan(PI*(Sum_E - .5));
        if((E_S0+resmass)>=minmass)
        {
            E = E_S0+resmass;
            
            if(resinfo->branchlist[0]->resinfo[0]->decay==true || resinfo->branchlist[0]->resinfo[1]->decay==true)
            {
                    double ma,mb,ma1,ma2,ma_pole,ma_0,ma_min,sum_ma,na,ma_gamma,ma_width;
                    double form_lambda,ma_kr,ma_k,ma_rho,ma_rho0,suma,rho_width,rho_width_0,spectsum,spectsum0,ma_kra,ma_ka,s0;
                
                    if(resinfo->branchlist[0]->resinfo[0]->decay==true)
                    {   ma_min=resinfo->branchlist[0]->resinfo[0]->minmass;
                        ma_pole=resinfo->branchlist[0]->resinfo[0]->mass;
                        mb=resinfo->branchlist[0]->resinfo[1]->mass;
                        ma_width=resinfo->branchlist[0]->resinfo[0]->width;
                        ma1=resinfo->branchlist[0]->resinfo[0]->branchlist[0]->resinfo[0]->mass;
                        ma2=resinfo->branchlist[0]->resinfo[0]->branchlist[0]->resinfo[1]->mass;
                        if(m1==776 && m2==138) { form_lambda=0.8; }
                        else if(resinfo->branchlist[0]->resinfo[1]->decay) { form_lambda=0.6; }
                        else if(resinfo->branchlist[0]->resinfo[0]->baryon==0) { form_lambda=1.6; }
                        else {form_lambda=2.0;}
                    }
                    if(resinfo->branchlist[0]->resinfo[1]->decay==true)
                    {   ma_min=resinfo->branchlist[0]->resinfo[1]->minmass;
                        ma_pole=resinfo->branchlist[0]->resinfo[1]->mass;
                        mb=resinfo->branchlist[0]->resinfo[0]->mass;
                        ma_width=resinfo->branchlist[0]->resinfo[1]->width;
                        ma1=resinfo->branchlist[0]->resinfo[1]->branchlist[0]->resinfo[0]->mass;
                        ma2=resinfo->branchlist[0]->resinfo[1]->branchlist[0]->resinfo[1]->mass;
                        if(m1==776 && m2==138) { form_lambda=0.8; }
                        else if(resinfo->branchlist[0]->resinfo[0]->decay) { form_lambda=0.6; }
                        else if(resinfo->branchlist[0]->resinfo[1]->baryon==0) { form_lambda=1.6; }
                        else {form_lambda=2.0;}
                    }
        
                    ma_kr=sqrt(abs(pow((ma_pole*ma_pole-ma1*ma1-ma2*ma2),2.0)-4.0*ma1*ma1*ma2*ma2))/(2.0*ma_pole);
                    suma=0.0;
                    int Na=100;
                
                    for(int na=0;na<Na;na++)
                    {
                        double sum_ma=(na+0.5)/Na;
                        ma_0 = 0.5*width*tan(PI*(sum_ma - .5));
                        ma = ma_0+ma_pole;
                        
                        if(ma>=ma_min && ma<=E)
                        {
                            ma_k=sqrt(abs(pow((ma*ma-ma1*ma1-ma2*ma2),2.0)-(4.0*ma1*ma1*ma2*ma2)))/(2.0*ma);
                            ma_gamma=ma_width*(ma_pole/ma)*((ma_k*ma_k*ma_k)/(ma_kr*ma_kr*ma_kr))*((ma_kr*ma_kr+HBARC*HBARC)/(ma_k*ma_k+HBARC*HBARC));
                            ma_rho=(2.0)/(ma_width*PI)*0.25*ma_gamma*ma_gamma/((0.25*ma_gamma*ma_gamma)+(ma_pole-ma)*(ma_pole-ma));
                            ma_rho0 = (1/PI)*(ma_width/2.0)/(0.25*ma_width*ma_width+ma_0*ma_0);
                            ma_kra=sqrt(abs(pow((resmass*resmass-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*resmass);
                            ma_ka=sqrt(abs(pow((E*E-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*E);
                            s0=ma+mb;
                            rho_width=(ma_ka*ma_ka*ma_ka)/(E*(ma_ka*ma_ka+HBARC*HBARC))*((pow(form_lambda,4.0)+0.25*pow((s0-resmass*resmass),2.0))/(pow(form_lambda,4.0)+pow((E*E-0.5*(s0+resmass*resmass)),2.0)));
                            rho_width_0=(ma_kra*ma_kra*ma_kra)/(resmass*(ma_kra*ma_kra+HBARC*HBARC));
                        
                            suma+=ma_rho/ma_rho0;
                            spectsum+=rho_width*ma_rho/ma_rho0;
                            spectsum0+=rho_width_0*ma_rho/ma_rho0;
                        }
                    }
                    if(ma>E) continue;
                    double avg_weight_ma=suma/Na;
                    double normal_ma=1.0/avg_weight_ma;
                    double spect=normal_ma*spectsum/Na;
                    double spect0=normal_ma*spectsum0/Na;
                    gamma=width*spect/spect0;
            }
            
            else{
                k=sqrt(abs(pow((E*E-m1*m1-m2*m2),2.0)-(4.0*m1*m1*m2*m2)))/(2.0*E);
                if(spin_deg<1.001)
                {gamma=width*(resmass/E)*(k/kr); }
                else {gamma=width*(resmass/E)*((k*k*k)/(kr*kr*kr))*((kr*kr+HBARC*HBARC)/(k*k+HBARC*HBARC));}
            }
            
            rho=(2.0)/(width*PI)*0.25*gamma*gamma/((0.25*gamma*gamma)+(resmass-E)*(resmass-E));
            rho_0 = (1/PI)*(width/2.0)/(0.25*width*width+E_S0*E_S0);
            freegascalc_onespecies(T,E,epsilon,P,dens,sigma2,dedt);
            weight=rho*dens/(rho_0*res_dens);
            if(weight>maxweight)
                maxweight=weight;
            esum+=epsilon*rho/rho_0;
            psum+=P*rho/rho_0;
            dsum+=dens*rho/rho_0;
            sigsum+=sigma2*rho/rho_0;
            dedtsum+=dedt*rho/rho_0;
            sum+=rho/rho_0;
        }
    }

    avg_weight=sum/N;
    normal=1.0/avg_weight;
    epsilon=normal*esum/N;
    P=normal*psum/N;
    dens=normal*dsum/N;
    sigma2=normal*sigsum/N;
    dedt=normal*dedtsum/N;
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
