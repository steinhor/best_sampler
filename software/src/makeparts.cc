#include "msu_sampler/sampler.h"
#include "msu_sampler/constants.h"

using namespace std;
using namespace msu_sampler;

int Csampler::MakeParts(Chyper *hyper){
	int nparts=0,dnparts,ires,nbose;
	CresInfo *resinfo;
	double udotdOmega=hyper->udotdOmega;
	double dN,dNtot=0,dNtotprime=0,mutot=0,I3;
	double dNcheck=0.0;
	double nptemp;
	CresMassMap::iterator iter;
	if(mastersampler->SETMU0)
		dNtot=dNtotprime=udotdOmega*nhadrons0;
	else
		dNtot=dNtotprime=udotdOmega*hyper->nhadrons;
	totvol+=udotdOmega;
	randy->netprob=0.0;
	randy->threshold=randy->ran_exp();
	if(randy->test_threshold(dNtot)){
		dNcheck=0.0;
		for(iter=reslist->massmap.begin();iter!=reslist->massmap.end();++iter){
			resinfo=iter->second;
			if(resinfo->pid!=22){
				nptemp=0;
				ires=resinfo->ires;
				I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
				mutot=hyper->muB*resinfo->baryon+hyper->muI*I3+hyper->muS*resinfo->strange;
				dN=exp(mutot)*density0i[ires]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				if(dNtotprime<-0.001){
					printf("dNtotprime=%g\n",dNtotprime);
					exit(1);
				}
				dnparts=CheckResInVolume(dN,Tf,resinfo,hyper);
				nparts+=dnparts;
				nptemp+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
		}
		if(bose_corr && !bose_test_off){
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				resinfo=reslist->GetResInfoPtr(211);
				nptemp=0;
				ires=resinfo->ires;
				mutot=nbose*hyper->muI;
				dN=exp(mutot)*pibose_dens0[nbose]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				nptemp+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
				resinfo=reslist->GetResInfoPtr(111);
				nptemp=0;
				ires=resinfo->ires;
				mutot=0.0;
				dN=pibose_dens0[nbose]*udotdOmega;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				nptemp+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
				resinfo=reslist->GetResInfoPtr(-211);
				nptemp=0;
				ires=resinfo->ires;
				mutot=-nbose*hyper->muI;
				dN=exp(mutot)*pibose_dens0[nbose]*udotdOmega;
				dNcheck+=dN;
				dNtotprime-=dN;
				dnparts=CheckResInVolume(dN,Tf/double(nbose),resinfo,hyper);
				nparts+=dnparts;
				nptemp+=dnparts;
				if(!(randy->test_threshold(dNtotprime))){
					randy->increment_netprob(dNtotprime);
					goto NoMoreParts;
				}
			}
		}
		NoMoreParts:
		if(dNcheck>dNtot*1.01){
			printf("Inside Csampler::MakeParts dNcheck=%g > dNtot=%g, dNtotprime=%g, T=%g\n",
			dNcheck,dNtot,dNtotprime,Tf);
			exit(1);
		}
		return nparts;
	}
	else{
		randy->increment_netprob(dNtot);
		return 0;
	}
}

int Csampler::CheckResInVolume(double dN,double T,CresInfo *resinfo,Chyper *hyper){
	int dnparts=0,alpha,ibin;
	double pmag,eta;
	FourVector p,r,ubj,ptilde,rtilde;
	randy->increment_netprob(dN);
	while(randy->test_threshold(0.0)){
		GetP(hyper,resinfo,p,T);
		if(BJORKEN_2D){
			if(fabs(hyper->u[3])>0.01){
				printf("BJORKEN_2D set, but u_z=%g\n",hyper->u[3]);
				exit(1);
			}
			eta=BJORKEN_YMAX*(1.0-randy->ran());
			ubj[1]=ubj[2]=0.0;
			ubj[0]=cosh(eta);
			ubj[3]=sinh(eta);
			Misc::Boost(ubj,p,ptilde);
			Misc::Boost(ubj,hyper->r,rtilde);
			for(alpha=0;alpha<4;alpha++){
				p[alpha]=ptilde[alpha];
				r[alpha]=rtilde[alpha];
			}
		}
		if (bose_test==true && (resinfo->pid==211||resinfo->pid==-211||resinfo->pid==111)) { //only runs when testing bose
      pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
      ibin=floorl(pmag/dp);
      dN_dp_p2[ibin]+=(2*M_PI*M_PI*HBARC*HBARC*HBARC)/(3*(ibin*dp+dp/2)*(ibin*dp+dp/2)*dp);
    }
		for(alpha=0;alpha<4;alpha++)
			r[alpha]=hyper->r[alpha];
		mastersampler->partlist->AddPart(resinfo->pid,p,r);
		dnparts+=1;
		randy->increase_threshold();
	}	return dnparts;
}

void Csampler::GetP(Chyper *hyper,CresInfo *resinfo,FourVector &p,double T){
	bool reflect;
	double pdotdOmega,nhatnorm,nhatdotp,wreflect;
	FourVector dOmegaTilde,ptilde;
	unsigned int alpha,beta;
	for (alpha=0;alpha<4;alpha++) {
		ptilde[alpha]=0;
	}
	FourVector pnoviscous;
	double m,nhat[4]={0.0};
	if((!resinfo->decay) || resinfo->width<0.0001 || USE_POLE_MASS){
		m=resinfo->mass;
		randy->generate_boltzmann(m,T,pnoviscous);
	}
	else{
		m=GenerateThermalMass(resinfo);
		randy->generate_boltzmann(m,T,pnoviscous);
	}
        printf("Likely a bug (undefined behaviour): viscous corrections = %d\n", VISCOUSCORRECTIONS);
	if(VISCOUSCORRECTIONS){
		for(alpha=1;alpha<4;alpha++){
			ptilde[alpha]=pnoviscous[alpha];
			for(beta=1;beta<4;beta++){
				if(mastersampler->SETMU0)
					ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((epsilon0+P0)*hyper->lambda*(M_PI));
				else
					ptilde[alpha]+=hyper->pitilde[alpha][beta]*pnoviscous[beta]/((hyper->epsilon+hyper->P)*hyper->lambda*(M_PI));
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
