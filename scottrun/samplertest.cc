#include "msu_sampler/master.h"
using namespace std;
using namespace msu_sampler;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid
int main(){
	const double PI=4.0*atan(1.0);
	vector<double> spectra_rho,spectra_delta;
	double dp=0.01,d3p,pmag;
	vector<double> massdist_rho,massdist_delta;
	double dm=0.005,norm;
	spectra_rho.resize(200);
	spectra_delta.resize(200);
	massdist_rho.resize(400);
	massdist_delta.resize(400);
	long long int npartstot=0;
	int nparts;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.txt");
	CpartList *partlist=new CpartList(&parmap);
	CmasterSampler ms(&parmap);
	ms.partlist=partlist;
	ms.randy->reset(time(NULL));

	// This makes a dummy hyper element for testing
	ms.MakeDummyHyper();
	Chyper *hyper=*(ms.hyperlist.begin());
	double V0=100000.0;
	hyper->T=0.145;
	hyper->T=0.5*parmap.getD("SAMPLER_TFMIN",0.150)+0.5*parmap.getD("SAMPLER_TFMAX",0.150);
	hyper->sigma=0.093;
	hyper->rhoB=0.0;
	hyper->rhoS=0.0;
	hyper->rhoI=0.0;
	hyper->epsilon=0.3;
	hyper->muB=hyper->muS=hyper->muI=0.0;
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
	// Now test particle production for specific pid
	long long unsigned int count=0;
	int pid;
	double totalenergy=0.0;
	printf("Enter pid to check: ");
	scanf("%d",&pid);
	CresInfo *resinfo=ms.reslist->GetResInfoPtr(pid);
	for(int ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		count+=partlist->CountResonances(pid);
		//partlist->IncrementSpectra(pid,dp,spectra);
		partlist->IncrementSpectra(113,dp,spectra_rho);
		partlist->IncrementSpectra(213,dp,spectra_rho);
		
		partlist->IncrementMassDist(113,dm,massdist_rho);
		partlist->IncrementMassDist(213,dm,massdist_rho);
		
		partlist->IncrementSpectra(1114,dp,spectra_delta);
		partlist->IncrementSpectra(2114,dp,spectra_delta);
		partlist->IncrementSpectra(2214,dp,spectra_delta);
		partlist->IncrementSpectra(2224,dp,spectra_delta);
		
		partlist->IncrementMassDist(1114,dm,massdist_delta);
		partlist->IncrementMassDist(2114,dm,massdist_delta);
		partlist->IncrementMassDist(2214,dm,massdist_delta);
		partlist->IncrementMassDist(2224,dm,massdist_delta);
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
	double ymax=parmap.getD("SAMPLER_BJORKEN_YMAX",0.001);
	printf("E/N=%g =? %g\n",EoverN,epsilonovern*sinh(ymax)/ymax);
	FILE *fptr;
	
	// First do rho
	if(parmap.getB("SAMPLER_USE_POLE_MASS",false))
		fptr=fopen("figs/rhospectra_polemass.dat","w");
	else
		fptr=fopen("figs/rhospectra.dat","w");
	for(unsigned int ip=0;ip<spectra_rho.size();ip++){
		pmag=(ip+0.5)*dp;
		d3p=4.0*PI*pmag*pmag*dp*hyper->udotdOmega*ms.NEVENTS_TOT;
		fprintf(fptr,"%6.3f %g\n",(ip+0.5)*dp,(1.0/3.0)*spectra_rho[ip]/d3p);
	}
	fclose(fptr);
	
	norm=0.0;
	for(unsigned int im=0;im<massdist_rho.size();im++)
		norm+=massdist_rho[im];
	if(parmap.getB("SAMPLER_USE_POLE_MASS",false))
		fptr=fopen("figs/rhomassdist_polemass.dat","w");
	else
		fptr=fopen("figs/rhomassdist.dat","w");
	for(unsigned int im=0;im<massdist_rho.size();im++){
		fprintf(fptr,"%6.3f %g\n",(im+0.5)*dm,massdist_rho[im]/(dm*norm));
	}
	fclose(fptr);
	
	// Now do Delta
	if(parmap.getB("SAMPLER_USE_POLE_MASS",false))
		fptr=fopen("figs/deltaspectra_polemass.dat","w");
	else
		fptr=fopen("figs/deltaspectra.dat","w");
	for(unsigned int ip=0;ip<spectra_delta.size();ip++){
		pmag=(ip+0.5)*dp;
		d3p=4.0*PI*pmag*pmag*dp*hyper->udotdOmega*ms.NEVENTS_TOT;
		fprintf(fptr,"%6.3f %g\n",(ip+0.5)*dp,0.125*spectra_delta[ip]/d3p);
	}
	fclose(fptr);
	
	norm=0.0;
	for(unsigned int im=0;im<massdist_delta.size();im++)
		norm+=massdist_delta[im];
	if(parmap.getB("SAMPLER_USE_POLE_MASS",false))
		fptr=fopen("figs/deltamassdist_polemass.dat","w");
	else
		fptr=fopen("figs/deltamassdist.dat","w");
	for(unsigned int im=0;im<massdist_delta.size();im++){
		fprintf(fptr,"%6.3f %g\n",(im+0.5)*dm,massdist_delta[im]/(dm*norm));
	}
	fclose(fptr);
	return 0;
}
