#ifndef __PART_CC__
#define __PART_CC__

#include "pratt_sampler/part.h"

CresList *Cpart::reslist=NULL;

Cpart::Cpart(){
	tau0=0.0;
	msquared=0.0;
	resinfo=NULL;
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
	y=eta=0.0;
}

Cpart::~Cpart(){
}

void Cpart::Copy(Cpart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	y=part->y;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		r[alpha]=part->r[alpha];
	}
	msquared=part->msquared;
	resinfo=part->resinfo;
}

void Cpart::CopyPositionInfo(Cpart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=part->r[alpha];
	}
}

void Cpart::CopyMomentumInfo(Cpart *part){  //copies all info except actionmap
	int alpha;
	y=part->y;
	msquared=part->msquared;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
	}
}

void Cpart::Init(int pid,double rxset,double ryset,double tauset,double etaset,double pxset,double pyset,double mset,double rapidityset){
	double et;
	resinfo=reslist->GetResInfoPtr(pid);
	p[1]=pxset; p[2]=pyset; msquared=mset*mset; y=rapidityset;
	r[1]=rxset; r[2]=ryset; tau0=tauset; eta=etaset;
	r[3]=tau0*sinh(eta);
	r[0]=tau0*cosh(eta);
	et=sqrt(p[1]*p[1]+p[2]*p[2]+msquared);
	p[3]=et*sinh(y);
	Setp0();
}

void Cpart::Print(){
	printf("________________ PART INFO FOR PART, pid=%d _____________________________\n",resinfo->code);
	printf("y=%g, Minv^2=%g,p=(%g,%g,%g,%g)\n",y,msquared,p[0],p[1],p[2],p[3]);
	printf("tau0=%g, eta=%g, r=(%g,%g,%g,%g)\n",tau0,eta,r[0],r[1],r[2],r[3]);
	printf("p=(%15.9e,%15.9e,%15.9e,%15.9e), y=%g =? %g\n",p[0],p[1],p[2],p[3],y,atanh(p[3]/p[0]));

	printf("________________________________________________________________________\n");
}

double Cpart::GetMass(){
	if(resinfo->code==22)
		return 0.0;
	else
		return sqrt(msquared);
}

void Cpart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
}

void Cpart::SetY(){
	y=asinh(p[3]/GetMT());
}

void Cpart::SetMass(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

double Cpart::GetEta(double tau){
	double dy,deta,dtau0,dtau;
	dy=y;
	deta=eta;
	dtau0=tau0;
	dtau=tau;
	deta=dy-asinh((dtau0/dtau)*sinh(dy-deta));
	return deta;
}

double Cpart::GetMT(){
	if(p[0]<fabs(p[3])){
		printf("Cpart::GetMT, catastrophe\n");
		Print();
		exit(1);
	}
	return sqrt(p[0]*p[0]-p[3]*p[3]);
}

double Cpart::GetPseudoRapidity(){
	double pmag,eta_ps;
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	eta_ps=atanh(p[3]/pmag);
	return eta_ps;
}

double Cpart::GetRapidity(){
	return 0.5*log((p[0]+p[3])/(p[0]-p[3]));
}

void Cpart::Boost(FourVector &u){
	BoostP(u);
	BoostR(u);
}

void Cpart::BoostP(FourVector &u){
	int alpha;
	FourVector pprime;
	Misc::Boost(u,p,pprime);
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pprime[alpha];
	y=atanh(p[3]/p[0]);
}

void Cpart::BoostR(FourVector &u){
	int alpha;
	FourVector rprime;
	Misc::Boost(u,r,rprime);
	for(alpha=0;alpha<4;alpha++)
		r[alpha]=rprime[alpha];
	eta=atanh(r[3]/r[0]);
	tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
}

void Cpart::GetHBTPars(double &t,double &rout,double &rside,double &rlong){
	const double tcompare=15.0;
	double pt,ptsquared,et;
	rlong=tau0*sinh(eta-y);
	t=tau0*cosh(eta-y);
	ptsquared=p[1]*p[1]+p[2]*p[2];
	pt=sqrt(ptsquared);
	et=sqrt(ptsquared+msquared);
	rout=(p[1]*r[1]+p[2]*r[2])/pt;
	rout=rout-(pt/et)*(t-tcompare);
	rside=(p[1]*r[2]-p[2]*r[1])/pt;
}

#endif
