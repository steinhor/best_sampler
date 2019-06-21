#include "sampler.h"
using namespace std;

int Csampler::MakeParts(Chyper *hyper){
	int nparts=0,ires;
	FourVector p;
	CresInfo *resinfo;
	double udotdOmega=hyper->udotdOmega; //Misc::DotProduct(hyper->u,hyper->dOmega);
	double dN,dNtot,dNprime,xx,mutot,I3;
	CresMassMap::iterator iter;
	if(mastersampler->SETMU0)
		dNtot=dNprime=udotdOmega*nhadronsf0;
	else
		dNtot=dNprime=udotdOmega*nhadronsf;
    //printf("randy->threshold=%g, randy->netprob=%g, dNtot=%g\n",randy->threshold,randy->netprob,dNtot);
	if(randy->test_threshold(dNtot)){
		ires=0;
		for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
			resinfo=iter->second;
			if(resinfo->code!=22){
				if(mastersampler->SETMU0==true)
					xx=1.0;
				else{
					I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
					mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
					xx=exp(mutot);
				}
				dN=densityf0[ires]*udotdOmega*xx;
				randy->increment_netprob(dN);
				dNprime-=dN;
				while(randy->test_threshold(0.0)){
					GetP(hyper,resinfo,p);
					nparts+=1;
					randy->increase_threshold();
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
		return nparts;
	}
	else{
		randy->increment_netprob(dNtot);
		return 0;
	}
}

void Csampler::GetP(Chyper *hyper,CresInfo *resinfo,FourVector &p){
	bool VISCOUSCORRECTIONS=true;
	double weight;
	bool success=false,reflect;
	int nparts=0;
	double pdotdOmega,eta,y,nhatnorm,nhatdotp,wreflect;
	FourVector dOmegaTilde,ptilde;
	int ispecies,alpha,beta,intweight,n,nsample;
	FourVector pnoviscous,u;
	double m,delN,r[3],w[3],nhat[4]={0.0};
	double mw;
	mw=maxweight[resinfo->ires];
	if(mw<0.0 || resinfo->width<0.001){
		m=resinfo->mass;
		randy->generate_boltzmann(m,Tf,pnoviscous);
	}
	else{
		m=GenerateThermalMass(resinfo);
		randy->generate_boltzmann(m,Tf,pnoviscous);
	}
	if(VISCOUSCORRECTIONS){
		for(alpha=1;alpha<4;alpha++){
			ptilde[alpha]=pnoviscous[alpha];
			for(beta=1;beta<4;beta++){
				if(mastersampler->SETMU0)
					ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((epsilonf0+Pf0)*lambdaf);
				else
					ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((epsilonf+Pf)*lambdaf);
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
	Misc::Boost(u,ptilde,p);
}
