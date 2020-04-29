#include "msu_sampler/master.h"
using namespace std;
using namespace msu_sampler;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid
int main(){
	long long int npartstot=0,nucleoncount=0,pioncount=0;
	double nucleonE=0.0,pionE=0.0;
	int nparts,alpha,beta,ievent;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.txt");
	CpartList *partlist=new CpartList(&parmap);
	CmasterSampler ms(&parmap);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));

	// This makes a dummy hyper element for testing
	ms.MakeDummyHyper();
	Chyper *hyper=*(ms.hyperlist.begin());
	double V0=10000.0,Tij;
	hyper->T=0.15;
	hyper->T=0.5*parmap.getD("SAMPLER_TFMIN",0.150)+0.5*parmap.getD("SAMPLER_TFMAX",0.150);
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
	
	printf("Enter pizz: ");
	scanf("%lf",&(hyper->pitilde[3][3]));
	hyper->pitilde[1][1]=hyper->pitilde[2][2]=-0.5*hyper->pitilde[3][3];
	hyper->CalcBiggestpitilde();
	
	(ms.ChooseSampler(hyper))->GetNHMu0();
	(ms.ChooseSampler(hyper))->GetMuNH(hyper);
	(ms.ChooseSampler(hyper))->CalcNHadronsEpsilonP(hyper);
	hyper->Print();
	printf("----- TARGET SE TENSOR /P -------\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			Tij=hyper->pitilde[alpha][beta];
			printf("%12.9f ",Tij/hyper->P);
		}
		printf("\n");
	}
	
	//hyper->muB=0.0502531705905960349327493288057178011854/hyper->T;  //exp(2*muB/T)=2
	
	for(ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		//partlist->IncrementSpectra(pid,dp,spectra);
		partlist->SumSETensor();
		nucleoncount+=partlist->CountResonances(-2212)+partlist->CountResonances(2212)+partlist->CountResonances(2112)+partlist->CountResonances(2112);
		pioncount+=partlist->CountResonances(211)+partlist->CountResonances(-211)+partlist->CountResonances(111);
		nucleonE+=partlist->SumEnergy(-2212)+partlist->SumEnergy(2212)+partlist->SumEnergy(2112)+partlist->SumEnergy(2112);
		pionE+=partlist->SumEnergy(211)+partlist->SumEnergy(-211)+partlist->SumEnergy(111);
		if((10*(ievent+1))%ms.NEVENTS_TOT==0){
			printf("finished %d percent\n",(ievent+1)*100/ms.NEVENTS_TOT);
		}
	}
	printf("rhoB=%g, rhoI=%g, rhoS=%g\n",hyper->rhoB,hyper->rhoI,hyper->rhoS);
	nucleonE=nucleonE/double(nucleoncount);
	pionE=pionE/double(pioncount);
	printf("npartstot=%lld =? %g, ratio=%g\n<E> for pions=%g, <E> for nucleons=%g\n",
	npartstot,hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT),
	double(npartstot)/(hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT)),pionE,nucleonE);
	printf("bulk momentum scaling ratio=%g\n",1.0/(1.0+hyper->PItilde/(3.0*hyper->Rbulk)));
	printf("------SE TENSOR -------\n");
	double norm=1.0/(hyper->udotdOmega*ms.NEVENTS_TOT);
	double PIbulk=norm*(partlist->SE[1][1]+partlist->SE[2][2]+partlist->SE[3][3])/3.0;
	PIbulk=PIbulk-hyper->P;
	//PIbulk=PIbulk+(partlist->SE[0][0]*norm-hyper->epsilon)*(hyper->P+hyper->epsilon)/(hyper->T*hyper->dedT);
	printf("T[0][0]=%g =? %g, ratio=%g\n",partlist->SE[0][0]*norm,hyper->epsilon,partlist->SE[0][0]*norm/hyper->epsilon);
	printf("PIbulk/P=%g =? %g, ratio=%g, P=%g, Tbulk=%g, deldotv*deltau=%g\n",
	PIbulk/hyper->P,hyper->PItilde/hyper->P,PIbulk/hyper->PItilde,hyper->P,hyper->T-hyper->PItilde/hyper->RTbulk,hyper->PItilde/hyper->Rbulk);
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			Tij=partlist->SE[alpha][beta]*norm;
			if(alpha==beta)
				Tij=Tij-hyper->P-PIbulk;
			printf("%12.9f ",Tij/hyper->P);
		}
		printf("\n");
	}
	
	char filename[80];
	sprintf(filename,"results/sheartest_rhoB%g.dat",hyper->rhoB);
	FILE *fptr=fopen(filename,"a");
	fprintf(fptr,"#pizz_hydro    pizz_h/P    pizz/pizz_h     e/etarget  \n");
	printf("%g, %g, ratio=%g=?1\n",(partlist->SE[3][3]*norm-hyper->P-PIbulk),hyper->pitilde[3][3],(partlist->SE[3][3]*norm-hyper->P-PIbulk)/hyper->pitilde[3][3]);
	fprintf(fptr,"%13.6e %13.6e %13.6e %13.6e\n",
	hyper->pitilde[3][3],hyper->pitilde[3][3]/hyper->P,(partlist->SE[3][3]*norm-hyper->P-PIbulk)/hyper->pitilde[3][3],partlist->SE[0][0]*norm/hyper->epsilon);
	fclose(fptr);
	return 0;
}
