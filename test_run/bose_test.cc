#include "pratt_sampler/test.h"
#include "pratt_sampler/sampler.h"

void calc_bose(CmasterSampler *ms,Csampler *sampler,Chyper *hyper, FILE *file);

int main(int argc, char *argv[]){
	int nprint=100,iprint=0;

	CparameterMap parmap;
	parmap.ReadParsFromFile("../run/parameters.dat");

	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	CpartList pl=CpartList(&parmap);
	ms.partlist=&pl;
	ms.randy->reset(time(NULL));

	ms.ReadHyper2D();

  Csampler *sampler=new Csampler(.150,0.093);
  Chyper *hyper=new Chyper();
  sampler->bose_test=true;

	hyper->muB=0.0;
	hyper->sigma=.0930;
	hyper->rhoB=0.007;
	hyper->epsilon=.121;
	hyper->rhoI=0.1*hyper->rhoB;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;
  hyper->T=.150;
  hyper->u[0]=1;
  hyper->u[1]=0;
  hyper->u[2]=0;
  hyper->u[3]=0;

  FILE *file=fopen("corrections.txt","w");
  calc_bose(&ms,sampler,hyper,file);
  fclose(file);

  sampler->bose_test_off=true;

  file=fopen("no_corrections.txt","w");
  calc_bose(&ms,sampler,hyper,file);
  fclose(file);

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}

void calc_bose(CmasterSampler *ms,Csampler *sampler,Chyper *hyper, FILE *file) {
  int N=50000*1000;
  int nparts=0;
  double totvol=0.0;
  vector<double> dN_dp_p2;

	sampler->GetTfMuNH(hyper);
  bool bose_off=sampler->bose_test_off;
	sampler=new Csampler(hyper->T,hyper->sigma);
  sampler->bose_test=true;
  if (bose_off) sampler->bose_test_off=true;
	sampler->GetTfMuNH(hyper);

  sampler->GetNHadronsEpsilonPF(hyper);
  sampler->CalcDensitiesMu0();
  hyper->lambda=sampler->CalcLambdaF(hyper);
  sampler->mastersampler=ms;

	printf("T=%lf\n",sampler->Tf);
	nparts=0;

	double ndp=400,pbound=4.5;
	sampler->dp=pbound/ndp;
  double dp=sampler->dp;
  double Tf=.144535,m=.138;
	int ibin;
	Cpart *part;

	for(int i=0;i<ndp;i++) {
		sampler->dN_dp_p2.push_back(0.0);
	}
	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
	}
	for (int ibin=0;ibin<ndp;ibin++) {
		fprintf(file,"%lf %lf\n",ibin*dp+dp/2,sampler->dN_dp_p2[ibin]/sampler->totvol);
	}

	printf("nparts=%d totvol=%lf dens=%lf nhadrons0=%lf\n",nparts,sampler->totvol,nparts/sampler->totvol,sampler->nhadrons0);

}
