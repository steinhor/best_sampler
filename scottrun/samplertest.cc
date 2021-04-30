#include "msu_sampler/master.h"
using namespace std;
using namespace msu_sampler;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid
int main(){
	double I3;
	CresInfo *resinfo;
	long long int npartstot=0,iprint=0,nprint;
	int nparts;
	CparameterMap parmap;
	FILE *fptr;
	parmap.ReadParsFromFile("parameters.txt");
	bool USE_POLE_MASS=parmap.getB("SAMPLER_USE_POLE_MASS",false);
	CmasterSampler ms(&parmap);
	CpartList *partlist=new CpartList(&parmap,ms.reslist);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));

	// This makes a dummy hyper element for testing
	ms.MakeDummyHyper();
	Chyper *hyper=*(ms.hyperlist.begin());
	double V0=1.0;
	hyper->T0=0.150;
	hyper->sigma=0.093;
	hyper->rhoB=0.0;
	hyper->rhoS=0.0;
	hyper->rhoI=0.0;
	hyper->epsilon=0.3;
	hyper->muB=hyper->muS=hyper->muI=0.0;
	
	//hyper->muB=0.0502531705905960349327493288057178011854/hyper->T0;  //exp(2*muB/T)=2
	
	hyper->u[1]=hyper->u[2]=hyper->u[3]=0.0;
	//hyper->u[1]=sqrt(3.0);
	hyper->u[0]=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]+hyper->u[3]*hyper->u[3]);
	hyper->r[1]=hyper->r[2]=hyper->r[3]=0.0;
	hyper->r[0]=10.0;
	for(int alpha=0;alpha<4;alpha++)
		hyper->dOmega[alpha]=V0*hyper->u[alpha];
	hyper->udotdOmega=hyper->u[0]*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]
		-hyper->u[2]*hyper->dOmega[2]-hyper->u[3]*hyper->dOmega[3];
	for(int alpha=0;alpha<4;alpha++)
		for(int beta=0;beta<4;beta++)
			hyper->pitilde[alpha][beta]=0.0;

	nprint=ms.NEVENTS_TOT/10;
	printf("ms.NEVENTS_TOT=%lld, nprint=%lld\n",ms.NEVENTS_TOT,nprint);
	for(long long int ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		partlist->CountResonances();
		iprint+=1;
		if(iprint==nprint){
			printf("finished %lld percent\n",((ievent+1)*100)/ms.NEVENTS_TOT);
			iprint=0;
		}
	}
	printf("npartstot=%lld =? %g, ratio=%g\n",
	npartstot,hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT),
	double(npartstot)/(hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT)));

	double epsilon,P,dens,sigma2,dedt,mass,ntarget;
	CresMassMap::iterator iter;
	fptr=fopen("rescount.dat","w");
	for(iter=ms.reslist->massmap.begin();iter!=ms.reslist->massmap.end();++iter){
		mass=iter->first;
		resinfo=iter->second;
		if(resinfo->pid!=22){
			if(resinfo->decay && !USE_POLE_MASS)
				EOS::freegascalc_onespecies_finitewidth(resinfo,hyper->T0,epsilon,P,dens,sigma2,dedt);
			else
				EOS::freegascalc_onespecies(hyper->T0,resinfo->mass,epsilon,P,dens,sigma2,dedt);
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
			dens*=resinfo->degen*exp(hyper->muB*resinfo->baryon+hyper->muI*I3+hyper->muS*resinfo->strange);
			ntarget=dens*double(ms.NEVENTS_TOT)*hyper->udotdOmega;
			fprintf(fptr,"%8d %8.4f %10.6f %12.6e %lld \n",resinfo->pid,mass,double(resinfo->count)/ntarget,ntarget,resinfo->count);
		}
	}
	fclose(fptr);

	return 0;
}
