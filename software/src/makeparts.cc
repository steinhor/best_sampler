#include "sampler.h"
using namespace std;

int Csampler::MakeParts(Chyper *hyper){
    int nparts=0,ires;
    CresInfo *resinfo;
    Cpart *part;
    double udotdOmega=hyper->udotdOmega;
    double dN=0,dNtot=0,dNprime=0,xx=1.0,mutot=0,I3;
    //double nptemp;

    CresMassMap::iterator iter;
    if(mastersampler->SETMU0)
        dNtot=dNprime=udotdOmega*nhadronsf0;
    else
        dNtot=dNprime=udotdOmega*nhadronsf;

    totvol+=udotdOmega;
    if(randy->test_threshold(dNtot)){
        ires=0;
        for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
            resinfo=iter->second;
            if(resinfo->code!=22) {
                if((resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) && parmap->getB("BOSE_CORR",false)) {
                    //nptemp=0;
                    for (int i=1;i<=parmap->getI("N_BOSE_CORR",1);i++) {
                        dN=npidens[resinfo->code][i-1]*udotdOmega;
                        randy->increment_netprob(dN);
                        dNprime-=dN;
                        while (randy->test_threshold(0.0)){
                            part=new Cpart();
                            GetP(hyper,resinfo,part->p,part,Tf/i);
                            nparts++;
                            //nptemp++;
                            randy->increase_threshold();
                        }
                    }
                }
                else {
                    if(mastersampler->SETMU0) dN=densityf0[ires]*udotdOmega;
                    else dN=densityf[ires]*udotdOmega;
                    randy->increment_netprob(dN);
                    dNprime-=dN;
                    //nptemp=0;
                    while(randy->test_threshold(0.0)){
                        part=new Cpart();
                        GetP(hyper,resinfo,part->p,part,Tf);
                        nparts++;
                        //nptemp++;
                        randy->increase_threshold();
                    }
                }

                //if (DensityMap.count(resinfo->code)==0) DensityMap.insert(pair<int,double>(resinfo->code,nptemp));
                //else DensityMap[resinfo->code]+=nptemp;

                if(!(randy->test_threshold(dNprime))){
                    randy->increment_netprob(dNprime);
                    goto NoMoreParts;
                }
                ires++;
            }
        }
    NoMoreParts:
        return nparts;
    }
    else{
        randy->increment_netprob(dNtot);
        return 0;
    }
}

void Csampler::GetP(Chyper *hyper,CresInfo *resinfo,FourVector &p,Cpart *part, double T){
    bool VISCOUSCORRECTIONS=true;
    bool reflect;
    double weight,pdotdOmega,nhatnorm,nhatdotp,wreflect;
    FourVector dOmegaTilde,ptilde;
    for (int i=0;i<4;i++) {
        ptilde[i]=0;
    }
    int alpha,beta;
    FourVector pnoviscous,u;
    double m,nhat[4]={0.0};
    double mw;
    mw=maxweight[resinfo->ires];
    if(mw<0.0 || resinfo->width<0.001){
        m=resinfo->mass;
        randy->generate_boltzmann(m,T,pnoviscous);
    }
    else{
        m=GenerateThermalMass(resinfo);
        randy->generate_boltzmann(m,T,pnoviscous);
    }
    part->msquared=m*m;
    if(VISCOUSCORRECTIONS){
        for(alpha=1;alpha<4;alpha++){
            ptilde[alpha]=pnoviscous[alpha];
            for(beta=1;beta<4;beta++){
                if(mastersampler->SETMU0)
                    ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((epsilonf0+Pf0)*lambdaf*(2*PI));
                else
                    ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((epsilonf+Pf)*lambdaf*(2*PI));
            }
        }
        ptilde[0]=sqrt(ptilde[1]*ptilde[1]+ptilde[2]*ptilde[2]+ptilde[3]*ptilde[3]+m*m);
    }
    else{
        for(alpha=0;alpha<4;alpha++)
            ptilde[alpha]=pnoviscous[alpha];
    }

    Misc::BoostToCM(hyper->u,hyper->dOmega,dOmegaTilde);  //dOmegaTilde is dOmega in fluid (u=0) frame

    pdotdOmega=ptilde[0]*dOmegaTilde[0]-ptilde[1]*dOmegaTilde[1]-ptilde[2]*dOmegaTilde[2];
    wreflect=pdotdOmega/(ptilde[0]*dOmegaTilde[0]);
    reflect=false;
    if(wreflect<0.0)
        reflect=true;
    if(wreflect<1.0){
        if(wreflect<randy->ran())
            reflect=true;
    }
    if(reflect){
        nhatnorm=sqrt(dOmegaTilde[1]*dOmegaTilde[1]+dOmegaTilde[2]*dOmegaTilde[2]);
        nhat[1]=dOmegaTilde[1]/nhatnorm;
        nhat[2]=dOmegaTilde[2]/nhatnorm;
        nhatdotp=nhat[1]*p[1]+nhat[2]*p[2];
        p[1]-=2.0*nhat[1]*nhatdotp;
        p[2]-=2.0*nhat[2]*nhatdotp;
    }
    Misc::Boost(hyper->u,ptilde,p);
}
