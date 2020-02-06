#include "pratt_sampler/master.h"
#include <cstring>
using namespace std;
using namespace pratt_sampler;

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

  int nparts=0;
	int N=pow(10,8);
	vector<Cpart>::iterator it;
	double totvol=0.0;

  double p[4];
	double T[4][4];
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			T[i][j]=0.0;
		}
	}
	Cpart part;
	double dens;
	int code;
	double temp;
	double psquared;
	double P;
	double Tcalc=0;

  //dummy hyper element
  Chyper *hyper=new Chyper();
	hyper->tau=12;
	hyper->muB=0.0;
	hyper->muI=0.0;
	hyper->muS=0.0;
	hyper->T=.144;
	hyper->sigma=0.0930;
	hyper->epsilon=.121;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;
	for (int i=0;i<4;i++) hyper->u[i]=0;
	hyper->dOmega[0]=0.089;
	hyper->dOmega[1]=-0.005;
	hyper->dOmega[2]=0.0078;
	hyper->dOmega[3]=0.0;
	hyper->pitilde[0][0]=0.0;
	hyper->pitilde[1][1]=0.0;
	hyper->pitilde[2][2]=0.0;
	hyper->pitilde[1][2]=hyper->pitilde[2][1]=.001;
	hyper->pitilde[3][3]=0.0;
	hyper->pitilde[3][1]=hyper->pitilde[1][3]=0.0;
	hyper->pitilde[3][2]=hyper->pitilde[2][3]=0.0;
	hyper->pitilde[0][1]=hyper->pitilde[1][0]=0.0;
	hyper->pitilde[0][2]=hyper->pitilde[2][0]=0.0;
	hyper->pitilde[0][3]=hyper->pitilde[3][0]=0.0;

	Csampler *sampler=new Csampler(hyper->T,hyper->sigma);
  sampler->mastersampler=&ms;
  sampler->viscous_test=true;

	hyper->rhoB=.007;
	hyper->rhoI=0.1*hyper->rhoB;
	sampler->GetTfMuNH(hyper);
	printf("hyper->T=%lf\tsampler->Tf=%lf\n",hyper->T,sampler->Tf);
  sampler->GetNHadronsEpsilonPF(hyper);
  sampler->CalcDensitiesMu0();
	hyper->lambda=sampler->CalcLambdaF(hyper);

	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
	}

	//calculate P and T from particles produced
  vector<Cpart> partvec=sampler->mastersampler->partlist->partvec;
  for(int k=0;k<nparts;k++) { //for (it=partvec.begin();it!=partvec.end();it++) {//for each species
    part=partvec[k];
		part.p[0]=sqrt(part.p[1]*part.p[1]+part.p[2]*part.p[2]+part.p[3]*part.p[3]-part.msquared);
		psquared=part.p[1]*part.p[1]+part.p[2]*part.p[2]+part.p[3]*part.p[3];

		if (!isnan(part.p[0])) P+=hyper->nhadrons*(psquared/(3*part.p[0]))/(nparts);
		else printf("psquared=%lf part.p[0]=%lf nparts=%d\n",psquared,part.p[0],nparts);

		for (int i=0;i<4;i++) { //row
			for (int j=i;j<4;j++) { //column
				temp=hyper->nhadrons*(part.p[i]*part.p[j]/part.p[0])/nparts;
				if (!isnan(part.p[0])) T[i][j]+=temp;
			}
		}
  }

	//print results
  printf("\n");
	printf("shear tensor from particles=\t\t\tshear tensor from input=\n");
	for (int i=1;i<4;i++) {
		for (int j=1;j<4;j++) {
			if (i==j) T[i][j]-=P;
			if (T[i][j]==0) T[i][j]=T[j][i];
			printf("%15lf",T[i][j]);
		}
		printf("\t");
		for (int j=1;j<4;j++) {
			printf("%15lf",hyper->pitilde[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("nparts=%d\n",nparts);
	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
