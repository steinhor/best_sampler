#include "sampler.h"
using namespace std;

int Csampler::MakeParts(Chyper *hyper){
    int nparts=0,ires;
    CresInfo *resinfo;
    double udotdOmega=hyper->udotdOmega; //Misc::DotProduct(hyper->u,hyper->dOmega);
    double dN=0,dNtot=0,dNprime=0,xx=1.0,mutot=0,I3;
    double nptemp;
    Cpart *part;
    double sum=0;
    double ndp=100;
    double p,pbound=4.5;
    CresMassMap::iterator iter;
    //printf("nhadronsf0=%lf nhadronsf=%lf\n",nhadronsf0,nhadronsf);
    if(mastersampler->SETMU0)
        dNtot=dNprime=udotdOmega*nhadronsf0;
    else
        dNtot=dNprime=udotdOmega*nhadronsf;

    ires=0;
    for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
        resinfo=iter->second;
        if(resinfo->code!=22){
            if(mastersampler->SETMU0) dN=densityf0[ires]*udotdOmega;
            else dN=densityf[ires]*udotdOmega;
            sum+=dN;
            //printf("dN=%lf densityf[ires]=%lf udotdOmega=%lf ires=%d xx=%lf resinfo->code=%d\n",dN,densityf[ires],udotdOmega,ires,xx,resinfo->code);
            ires++;
        }
    }
    //printf("dNtot=%lf\tsum=%lf\tmuB=%lf\tmuI=%lf\n",dNtot,sum,muB,muI);
    //if (hyper->ihyp==0) printf("%lf\n",sum-dNtot);

    //sum=0;
    totvol+=udotdOmega;
    if(randy->test_threshold(dNtot)){
        ires=0;
        for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
            resinfo=iter->second;
            if(resinfo->code!=22) {
                if(mastersampler->SETMU0) dN=densityf0[ires]*udotdOmega;
                else dN=densityf[ires]*udotdOmega;

                randy->increment_netprob(dN);
                dNprime-=dN;
                nptemp=0;
                while(randy->test_threshold(0.0)){
                    part=new Cpart();
                    GetP(hyper,resinfo,part->p,part);
                    nparts++;
                    nptemp++;
                    pmap[resinfo->code].push_back(part);
                    p=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3]);

                    if (resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) {
                    	for (double i=0;i<=pbound;i+=pbound/double(ndp)) {
                    		if (p>=i && p<(i+pbound/ndp)) {
                                if (dpmap.count(i+pbound/(2*ndp))==0) dpmap.insert(pair<double,double>(i+pbound/(2*ndp),(2*PI*PI*HBARC*HBARC*HBARC)/(p*p)));
                                else dpmap[i+pbound/(2*ndp)]+=1/(p*p);
                                break;
                            }
                    	}
                    }

                    randy->increase_threshold();
                }

                if (DensityMap.count(resinfo->code)==0) {
                    DensityMap.insert(pair<int,double>(resinfo->code,nptemp));
                }
                else {
                    DensityMap[resinfo->code]+=nptemp;
                }

                if(!(randy->test_threshold(dNprime))){
                    randy->increment_netprob(dNprime);
                    goto NoMoreParts;
                }
                ires++;
            }
            //printf("resinfo: ires=%d code=%d mass=%lf\n",resinfo->ires, resinfo->code, resinfo->mass);
        }
    NoMoreParts:
        //printf("dNtot=%lf\tsum=%lf\tdNprime=%lf\n",dNtot,sum,dNprime);
        return nparts;
    }
    else{
        randy->increment_netprob(dNtot);
        return 0;
    }
}

void Csampler::GetP(Chyper *hyper,CresInfo *resinfo,FourVector &p,Cpart *part){
    bool VISCOUSCORRECTIONS=true;
    double weight;
    bool success=false,reflect;
    int nparts=0;
    double pdotdOmega,eta,y,nhatnorm,nhatdotp,wreflect;
    double temp;
    FourVector dOmegaTilde,ptilde;
    int ispecies,alpha,beta,intweight,n,nsample;
    FourVector pnoviscous,u;
    double m,delN,r[3],w[3],nhat[4]={0.0};
    double mw;
    mw=maxweight[resinfo->ires];
    //printf("test 9\n");
    if(mw<0.0 || resinfo->width<0.001){
        m=resinfo->mass;
        randy->generate_boltzmann(m,Tf,pnoviscous);
    }
    else{
        m=GenerateThermalMass(resinfo);
        randy->generate_boltzmann(m,Tf,pnoviscous);
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
    //printf("ptilde=(%g,%g,%g,%g)\n",ptilde[0],ptilde[1],ptilde[2],ptilde[3]);

    Misc::BoostToCM(hyper->u,hyper->dOmega,dOmegaTilde);  //dOmegaTilde is dOmega in fluid (u=0) frame
    pdotdOmega=ptilde[0]*dOmegaTilde[0]-ptilde[1]*dOmegaTilde[1]-ptilde[2]*dOmegaTilde[2];
    //printf("u=(%g,%g), dOmega=(%g,%g,%g), dOmegaTilde=(%g, %g,%g)\n",
    //ux,uy,dOmega[0],dOmega[1],dOmega[2],dOmegaTilde[0],dOmegaTilde[1],dOmegaTilde[2]);
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
    for (int i=0;i<4;i++) {
        //printf("ptilde=%lf p[%d]=%lf\t",ptilde[i],i,p[i]);
    }
    Misc::Boost(hyper->u,ptilde,p);
    for (int i=0;i<4;i++) {
        //printf("ptilde=%lf p[%d]=%lf\t",ptilde[i],i,p[i]);
    }
    //printf("\n\n");
}
