#include "msu_sampler/master.h"
using namespace std;
using namespace msu_sampler;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid
int main(){
	double T0=0.150,V0=2.0;
	long long int npartstot=0,iprint=0,nprint,ievent;
	int nparts,nhyper=100000;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.txt");
	CmasterSampler ms(&parmap);
	CpartList *partlist=new CpartList(&parmap,ms.reslist);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));
	nprint=ms.NEVENTS_TOT/10;

	// This makes a dummy hyper element for testing
	Chyper *hyper;
	list<Chyper *>::iterator it_hyper;
		
	// Make some dummy hyper elements
	ms.MakeDummyHyper(nhyper);
	for(ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		for(it_hyper=ms.hyperlist.begin();it_hyper!=ms.hyperlist.end();it_hyper++){
			hyper=*it_hyper;
			hyper->T0=T0;
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
				hyper->dOmega[alpha]=V0*hyper->u[alpha]*2.0*ms.randy->ran();
			hyper->udotdOmega=hyper->u[0]*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]
				-hyper->u[2]*hyper->dOmega[2]-hyper->u[3]*hyper->dOmega[3];
			for(int alpha=0;alpha<4;alpha++)
				for(int beta=0;beta<4;beta++)
					hyper->pitilde[alpha][beta]=0.0;
		}
		// Perform events
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		iprint+=1;
		if(iprint==nprint){
			printf("finished %lld percent\n",((ievent+1)*100)/ms.NEVENTS_TOT);
			iprint=0;
			it_hyper=ms.hyperlist.begin();
			hyper=*it_hyper;
			printf("npartstot=%lld =? %g, ratio=%g\n",
			npartstot,nhyper*hyper->nhadrons*V0*double(ms.NEVENTS_TOT),
			double(npartstot)/(nhyper*hyper->nhadrons*V0*double(ms.NEVENTS_TOT)));
		}
	}
	ms.ClearHyperList();
	
	delete partlist;
	return 0;
}
