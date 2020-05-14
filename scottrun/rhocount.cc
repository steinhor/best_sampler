#include "msu_sampler/master.h"
#include <cstring>
using namespace std;
using namespace msu_sampler;
int main(){
	double T,V0=1000.0;
	long long int deltacount=0,rhocount=0,npartstot=0;
	int nparts,alpha,beta,ievent,nevents;
	CparameterMap parmap;
	parmap.ReadParsFromFile("hydroparameters.txt");
	CpartList *partlist=new CpartList(&parmap);
	CmasterSampler ms(&parmap);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));
	FILE *fptr=fopen("rhocount.dat","w");

	// This makes a dummy hyper element for testing
	ms.MakeDummyHyper();
	Chyper *hyper=*(ms.hyperlist.begin());
	hyper->sigma=0.093;
	hyper->rhoB=0.0;
	hyper->rhoS=0.0;
	hyper->rhoI=0.0;
	hyper->muB=0.0;
	hyper->muS=hyper->muI=0.0;
	hyper->u[1]=hyper->u[2]=hyper->u[3]=0.0;
	//hyper->u[1]=sqrt(3.0);
	hyper->u[0]=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]+hyper->u[3]*hyper->u[3]);
	hyper->r[1]=hyper->r[2]=hyper->r[3]=0.0;
	hyper->r[0]=10.0;
	for(alpha=0;alpha<4;alpha++)
		hyper->dOmega[alpha]=V0*hyper->u[alpha];
	hyper->udotdOmega=hyper->u[0]*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]
		-hyper->u[2]*hyper->dOmega[2]-hyper->u[3]*hyper->dOmega[3];
	for(alpha=0;alpha<4;alpha++)
		for(beta=0;beta<4;beta++)
			hyper->pitilde[alpha][beta]=0.0;
	hyper->PItilde=0.0;
	nevents=ms.NEVENTS_TOT;
	
	for(T=100;T<161;T+=5.0){
		printf("------- T=%g ---------\n",T);
		hyper->T=T/1000.0;
		npartstot=deltacount=rhocount=0;
		hyper->CalcBiggestpitilde();
		(ms.ChooseSampler(hyper))->GetNHMu0();
		(ms.ChooseSampler(hyper))->GetMuNH(hyper);
		(ms.ChooseSampler(hyper))->CalcNHadronsEpsilonP(hyper);
		printf("Tf=%g\n",(ms.ChooseSampler(hyper))->Tf);
	
		for(ievent=0;ievent<nevents;ievent++){
			nparts=ms.MakeEvent();
			npartstot+=nparts;
			deltacount+=ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
				+ms.partlist->CountResonances(2214)+ms.partlist->CountResonances(-2214)
					+ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
						+ms.partlist->CountResonances(1114)+ms.partlist->CountResonances(-1114);
			rhocount+=ms.partlist->CountResonances(113)+ms.partlist->CountResonances(213)+ms.partlist->CountResonances(-213);
			if((10*(ievent+1))%nevents==0)
				printf("ievent=%d\n",ms.NEVENTS);
		}
		printf("nparts/event=%g\n",double(npartstot)/double(nevents));

		printf("rho_rho=%g\n",double(rhocount)/(V0*nevents));
		printf("rho_delta=%g\n",double(deltacount)/(V0*nevents));
		fprintf(fptr,"%6.1f %10.5f %10.5f\n",1000*hyper->T,double(rhocount)/(V0*nevents),double(deltacount)/(V0*nevents));
	}
	fclose(fptr);
	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
