#include "pratt_sampler/master.h"
using namespace std;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid
int main(int argc, char *argv[]){
	long long int npartstot=0;
	int nparts;
	CparameterMap parmap;
	parmap.ReadParsFromFile("testparameters.txt");
	CpartList *partlist=new CpartList(&parmap);
	CmasterSampler ms(&parmap);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));
	
	// This makes a dummy hyper element for testing
	ms.MakeDummyHyper();
	Chyper *hyper=*(ms.hyperlist.begin());
	double V0=100.0;
	hyper->T=0.150;
	hyper->T=0.5*parmap.getD("SAMPLER_TFMIN",0.150)+0.5*parmap.getD("SAMPLER_TFMAX",0.150);
	hyper->sigma=0.093;
	hyper->rhoB=0.0;
	hyper->rhoS=0.0;
	hyper->rhoI=0.0;
	hyper->epsilon=0.3;
	hyper->muB=hyper->muS=hyper->muI=0.0;
	hyper->u[1]=hyper->u[2]=hyper->u[3]=0.0;
	hyper->u[1]=sqrt(3.0);
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
	// Now test particle production for specific pid
	long long int count=0;
	int pid;
	double totalenergy=0.0;
	printf("Enter pid to check: ");
	scanf("%d",&pid);
	CresInfo *resinfo=ms.reslist->GetResInfoPtr(pid);
	for(int ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		count+=partlist->CountResonances(pid);
		totalenergy+=partlist->SumEnergy(pid);
		if((10*(ievent+1))%ms.NEVENTS_TOT==0)
			printf("finished %d percent\n",(ievent+1)*100/ms.NEVENTS_TOT);
	}
	printf("npartstot=%lld =? %g, ratio=%g\n",
	npartstot,hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT),
	double(npartstot)/(hyper->nhadrons*hyper->udotdOmega*double(ms.NEVENTS_TOT)));
	
	double I3,epsilon,P,dens,sigma2,dedt,ntarget;
	if(resinfo->decay)
		EOS::freegascalc_onespecies_finitewidth(resinfo,hyper->T,epsilon,P,dens,sigma2,dedt);
	else
		EOS::freegascalc_onespecies(hyper->T,resinfo->mass,epsilon,P,dens,sigma2,dedt);
	I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
	dens*=resinfo->degen*exp(hyper->muB*resinfo->baryon+hyper->muI*I3+hyper->muS*resinfo->strange);
	ntarget=dens*double(ms.NEVENTS_TOT)*hyper->udotdOmega;
	printf("------ For PID=%d ------\n",pid);
	resinfo->PrintBranchInfo();
	printf("count=%lld =? %g, ratio=%g\n",count,ntarget,double(count)/ntarget);
	double EoverN=totalenergy/double(count);
	double epsilonovern=hyper->u[0]*epsilon*resinfo->degen/dens;
	printf("E/N=%g =? %g\n",EoverN,epsilonovern);
				
	return 0;
}