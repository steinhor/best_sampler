#include "master.h"
using namespace std;

int main(int argc, char *argv[]){
	int nprint=100,iprint=0;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	Csampler *sampler=new Csampler(0.150,0.0930);
	Crandy *randy=new Crandy(1234);
	ms.randy->reset(time(NULL));

	ms.ReadHyper2D();
	const int NHYPER=15;

	Chyper *hyper=new Chyper();
	Eigen::MatrixXd A(4,4);
	Eigen::MatrixXd B(4,4);
	Eigen::VectorXd dmu(4);
	dmu[0]=0.05; dmu[1]=dmu[2]=dmu[3]=0.002;
	Eigen::VectorXd mu(4);
	Eigen::VectorXd rho(4),rho1(4),rho2(4);
	double rhoB,rhoI,rhoS,epsilon;


	/* // For testing---
	double Tf=120.0;
	double muB=12.0,muI=0.0,muS=0.0;
	mu[0]=Tf; mu[1]=muB; mu[2]=muI; mu[2]=mu[3];
	sampler->Tf=mu[0]; sampler->muB=mu[1]; sampler->muI=mu[2]; sampler->muS=mu[3];
	sampler->GetEpsilonRhoDerivatives(epsilon,rhoB,rhoI,rhoS,A);
	rho[0]=epsilon; rho[1]=rhoB; rho[2]=rhoI; rho[3]=rhoS;
	printf("epsilon=%g, rhoB=%g, rhoI=%g, rhoS=%g\n",epsilon,rhoB,rhoI,rhoS);
	cout << A << endl;
	printf("-----------------------------\n");

	int i,j;
	for(j=0;j<4;j++){
		mu[j]+=0.5*dmu[j];
		sampler->Tf=mu[0]; sampler->muB=mu[1]; sampler->muI=mu[2]; sampler->muS=mu[3];
		sampler->GetEpsilonRhoDerivatives(epsilon,rhoB,rhoI,rhoS,A);
		rho2[0]=epsilon; rho2[1]=rhoB; rho2[2]=rhoI; rho2[3]=rhoS;
		mu[j]-=dmu[j];
		sampler->Tf=mu[0]; sampler->muB=mu[1]; sampler->muI=mu[2]; sampler->muS=mu[3];
		sampler->GetEpsilonRhoDerivatives(epsilon,rhoB,rhoI,rhoS,A);
		rho1[0]=epsilon; rho1[1]=rhoB; rho1[2]=rhoI; rho1[3]=rhoS;
		mu[j]+=0.5*dmu[j];
		for(i=0;i<4;i++)
			B(i,j)=(rho2[i]-rho1[i])/dmu[j];
	}
	cout << B << endl;
	printf("-----------------------------\n");
	//exit(1); */

	int nparts=0;
	int temp;
	list<Chyper *>::iterator it;
	map<int,double> MasterDensityMap;
	double totvol=0.0;
	for (int i=0; i<493; i++) { //initialize master map
		MasterDensityMap.insert(pair<int,double>(i,0.0));
	}
	for (int i=0;i<100;i++) { //for each sampler
		sampler=new Csampler(0.150,0.0930);
		sampler->totvol=0.0;
		for (int i=0; i<493; i++) { //initialize map
			sampler->DensityMap.insert(pair<int,double>(i,0.0));
		}
		for (it=ms.hyperlist.begin();it!=ms.hyperlist.end();it++) { //for each hyper element
			hyper=*it;
			temp=sampler->MakeParts(hyper);
			nparts+=temp;
		}

		//update master data
		totvol+=sampler->totvol;
		for (int i=0; i<493; i++) {
			MasterDensityMap[i]+=sampler->DensityMap[i];
		}
	}
	printf("ires:\tDensityMap:\tdensityf0:\n");
	for (int i=0; i<493; i++) {
		MasterDensityMap[i]=MasterDensityMap[i]/totvol;
		printf("%d\t%lf\t%lf\n",i,MasterDensityMap[i],sampler->densityf0[i]);
	}
	printf("nparts=%d\n",nparts);
	printf("totvol=%lf\n",totvol);
	printf("nhad=%lf",sampler->nhadronsf0);

	/*
	for(int ihyper=0;ihyper<=NHYPER;ihyper++){
		hyper->muB=0.0;
		hyper->T=80.0;
		hyper->sigma=93.0;
		for(int jhyper=0;jhyper<=NHYPER;jhyper++){
			hyper->rhoB=jhyper*0.15/double(NHYPER);
			hyper->epsilon=50+1000*hyper->rhoB+750.0*ihyper/double(NHYPER);
			hyper->rhoI=0.1*hyper->rhoB;
			hyper->rhoS=0.0;

			sampler->GetTfMuNH(hyper);
			printf("epsilon=%g, rhoB=%6.4f, nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n",hyper->epsilon,hyper->rhoB,hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
		}
		printf("--------------------------------------------------------------------------------------------\n");
	}
	*/


	printf("------ hyper elements created ----------- \n");

	/*
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		//printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	*/
	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
