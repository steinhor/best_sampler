#ifndef __PART_CC__
#define __PART_CC__

#include "msu_sampler/part.h"
using namespace msu_sampler;

Cpart::Cpart(){
	msquared=0.0;
	pid=0;
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
}

Cpart::~Cpart(){
}

void Cpart::Print(){
	printf("________________ PART INFO FOR PART, pid=%d _____________________________\n",pid);
	printf("m^2=%g\\p=(%g,%g,%g,%g)\n",msquared,p[0],p[1],p[2],p[3]);
	printf("r=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
}

double Cpart::GetMass(){
	if(pid==22)
		return 0.0;
	else
		return sqrt(msquared);
}

void Cpart::SetMsquared(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

void Cpart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
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
}

void Cpart::BoostR(FourVector &u){
	int alpha;
	FourVector rprime;
	Misc::Boost(u,r,rprime);
	for(alpha=0;alpha<4;alpha++)
		r[alpha]=rprime[alpha];
}

CpartList::CpartList(CparameterMap *parmap){
	nparts_blocksize=parmap->getI("SAMPLER_NPARTS_BLOCKSIZE",2000);
	partvec.resize(nparts_blocksize);
	Reset();
}

CpartList::~CpartList(){
	partvec.clear();
	partvec.shrink_to_fit();
}

Cpart* CpartList::GetPart(){
	if(partvec.size()==nparts){
		printf("resizing partvec, old size=%lu\n",partvec.size());
		partvec.resize(partvec.size()+nparts);
	}
	nparts+=1;
	return &partvec[nparts];
}
void CpartList::Kill(){
	partvec.clear();
	partvec.shrink_to_fit();
	nparts=0;
}

void CpartList::Reset(){
	nparts=0;
	for(int alpha=0;alpha<4;alpha++){
		for(int beta=0;beta<4;beta++)
			SE[alpha][beta]=0.0;
	}
}

void CpartList::WriteParts(string filename){
	FILE *fptr=fopen(filename.c_str(),"w");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open file to write parts\n");
		printf("Error %d \n", errno);
		exit(1);
	}
	for(unsigned int ipart=0;ipart<nparts;ipart++){
		fprintf(fptr,"%5d %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e\n",
		partvec[ipart].pid,partvec[ipart].msquared,partvec[ipart].p[0],partvec[ipart].p[1],partvec[ipart].p[2],partvec[ipart].p[3],
		partvec[ipart].r[0],partvec[ipart].r[1],partvec[ipart].r[2],partvec[ipart].r[3]);
	}
	fclose(fptr);
}

unsigned int CpartList::CountResonances(int pid){
	unsigned int count=0,ipart;
	for(ipart=0;ipart<nparts;ipart++){
		if(partvec[ipart].pid==pid)
			count+=1;
	}
	return count;
}

void CpartList::IncrementSpectra(int pid,double dp,vector<double> &spectra){
	unsigned int ip,ipart,np=spectra.size();
	double pt;
	for(ipart=0;ipart<nparts;ipart++){
		if(abs(partvec[ipart].pid)==pid){
			//pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]);
			pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]
				+partvec[ipart].p[3]*partvec[ipart].p[3]);
			ip=lrint(floor(pt/dp));
			if(ip<np)
				spectra[ip]+=1;
		}
	}
}

void CpartList::IncrementMassDist(int pid,double dm,vector<double> &massdist){
	unsigned int im,ipart,nm=massdist.size();
	double m;
	for(ipart=0;ipart<nparts;ipart++){
		if(abs(partvec[ipart].pid)==pid){
			partvec[ipart].SetMsquared();
			m=sqrt(partvec[ipart].msquared);
			im=lrint(floor(m/dm));
			if(im<nm)
				massdist[im]+=1;
		}
	}
}

double CpartList::SumEnergy(){
	double energy=0.0,ipart;
	for(ipart=0;ipart<nparts;ipart++){
		energy+=partvec[ipart].p[0];
	}
	return energy;
}

double CpartList::SumEnergy(int pid){
	double energy=0.0;
	unsigned int ipart;
	for(ipart=0;ipart<nparts;ipart++){
		if(partvec[ipart].pid==pid)
			energy+=partvec[ipart].p[0];
	}
	return energy;
}

void CpartList::AddPart(int pidset,FourVector &pset,FourVector &rset){
	unsigned int alpha;
	if(partvec.size()==nparts){
		//printf("resizing partvec, old size=%lu\n",partvec.size());
		partvec.resize(partvec.size()+nparts);
	}
	partvec[nparts].pid=pidset;
	for(alpha=0;alpha<4;alpha++){
		partvec[nparts].p[alpha]=pset[alpha];
		partvec[nparts].r[alpha]=rset[alpha];
	}
	partvec[nparts].SetMsquared();
	nparts+=1;
}

void CpartList::SumSETensor(){
	unsigned int ipart,alpha,beta;
	//nparts=partvec.size();
	for(ipart=0;ipart<nparts;ipart++){
		//printf("partvec[%d]=(%g,%g,%g,%g)\n",ipart,partvec[ipart].p[0],partvec[ipart].p[1],partvec[ipart].p[2],partvec[ipart].p[3]);
		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				SE[alpha][beta]+=partvec[ipart].p[alpha]*partvec[ipart].p[beta]/partvec[ipart].p[0];
			}
		}
	}
}

#endif
