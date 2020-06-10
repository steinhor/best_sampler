#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "msu_sampler/sampler.h"
#include "msu_sampler/constants.h"
#include "msu_sampler/misc.h"
#include <iostream>

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
				if(dNtotprime<-0.0001){
					printf("dNtotprime=%g, should not be negative, mutot=%g, dNcheck=%g, dNtot=%g\n",
					dNtotprime,mutot,dNcheck,dNtot);
					printf("nhadrons0=%g, hyper->nhadrons=%g\n",nhadrons0,hyper->nhadrons);
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
		GetP(hyper,T,resinfo,p);
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

void Csampler::GetPInFluidFrame(double m,Chyper *hyper,double T,FourVector &p){
	FourVector pnoshear,pnobulk;
	double Tbulk;
	int alpha;
	if(hyper->Rvisc_calculated==false){
		CalcRvisc(hyper);
		hyper->Rvisc_calculated=true;
	}
	if(INCLUDE_BULK_VISCOSITY){
		Tbulk=T-hyper->PItilde/hyper->RTbulk;
		randy->generate_boltzmann(m,Tbulk,pnobulk);
		//printf("pnobulk=(%g,%g,%g,%g)\n",pnobulk[0],pnobulk[1],pnobulk[2],pnobulk[3]);
		BulkScale(hyper,m,pnobulk,pnoshear);
		//for(alpha=0;alpha<4;alpha++)
		//pnoshear[alpha]=pnobulk[alpha];
	}
	else
		randy->generate_boltzmann(m,T,pnoshear);
	if(INCLUDE_SHEAR_VISCOSITY){
		ShearScale(hyper,m,pnoshear,p);
	}
	else{
		for(alpha=0;alpha<4;alpha++)
			p[alpha]=pnoshear[alpha];
	}
}

void Csampler::BulkScale(Chyper *hyper,double mass,FourVector &pnobulk,FourVector &p){
	int alpha;
	for(alpha=1;alpha<4;alpha++)
		p[alpha]=pnobulk[alpha]/(1.0+hyper->PItilde/(3.0*hyper->Rbulk));
	p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
}

void Csampler::ShearScale(Chyper *hyper,double mass,FourVector &pnoshear,FourVector &p){
	int alpha,beta;
	double R,Rmax,pmag,ctheta,stheta,phi;
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pnoshear[alpha];
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	Rmax=1.0-hyper->biggestpitilde*pmag*pmag/(p[0]*Tf*hyper->Rshear);
	//printf("pmag=%g, Rmax=%g, biggestpitilde=%g, Rshear=%g\n",pmag,Rmax,hyper->biggestpitilde,hyper->Rshear);

	R=1.0;
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			R-=hyper->pitilde[alpha][beta]*p[alpha]*p[beta]/(p[0]*Tf*hyper->Rshear);
		}
	}
	while(randy->ran()>R/Rmax){
		ctheta=1.0-2.0*randy->ran();
		stheta=sqrt(1.0-ctheta*ctheta);
		phi=2.0*M_PI*randy->ran();
		p[1]=pmag*stheta*cos(phi);
		p[2]=pmag*stheta*sin(phi);
		p[3]=pmag*ctheta;
		p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
		R=1.0;
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++){
				R-=hyper->pitilde[alpha][beta]*p[alpha]*p[beta]/(p[0]*Tf*hyper->Rshear);
			}
		}
		if(R>Rmax){
			printf("R=%g, pmag=%g, Rmax=%g, biggestpitilde=%g, Rshear=%g\n",R,pmag,Rmax,hyper->biggestpitilde,hyper->Rshear);
		}
	};
	
	/*
	for(alpha=1;alpha<4;alpha++){
		p[alpha]=pnoshear[alpha];
		for(beta=1;beta<4;beta++){
			p[alpha]-=hyper->pitilde[alpha][beta]*pnoshear[beta]/hyper->Rshear;
		}
	}
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
	*/
	
}

void Csampler::GetP(Chyper *hyper,double T,CresInfo *resinfo,FourVector &p){
	bool reflect;
	double pdotdOmega,nhatnorm,nhatdotp,wreflect;
	FourVector dOmegaTilde,ptilde;
	double m,nhat[4]={0.0};
	if((!resinfo->decay) || resinfo->width<0.0001 || USE_POLE_MASS)
		m=resinfo->mass;
	else
		m=GenerateThermalMass(resinfo);
	GetPInFluidFrame(m,hyper,T,ptilde);

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
