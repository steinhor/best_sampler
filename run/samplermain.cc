#include "master.h"
using namespace std;

int main(int argc, char *argv[]){
	int nprint=100,iprint=0;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	Csampler *sampler=new Csampler(0.150,0.0930);
	sampler->partmap=new CpartMap;
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

	//This is Scott's, I don't know what it's for.
	/*****************************************************************************************
	// For testing---
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
	//exit(1);
	****************************************************************************************/

	//Test whether the densities of particles produced in MakeParts match the expected densities.
	/****************************************************************************************/
	int nparts=0;
	int temp;
	int n=0;
	double avg_rhoI=0.0;
	double real_rhoI=0.0;
	list<Chyper *>::iterator it;
	map<int,double>::iterator itr;
	map<int,double> MasterDensityMap;
	double totvol=0.0;

	for (int i=0;i<100;i++) { //for each sampler
		if (i!=0) {
			sampler=new Csampler(0.150,0.0930);
			sampler->partmap=new CpartMap;
		}
		sampler->totvol=0.0;
		for (it=ms.hyperlist.begin();it!=ms.hyperlist.end();it++) { //for each hyper element
			hyper=*it;
			if (it==ms.hyperlist.begin()) sampler->GetTfMuNH(hyper);
			temp=sampler->MakeParts(hyper);
			nparts+=temp;
			n++;
		}
		//update master data
		totvol+=sampler->totvol;
		for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
			int code=itr->first;
			double dens=itr->second;
			if (MasterDensityMap.count(code)==0) {
				MasterDensityMap.insert(pair<int,double>(code,dens));
			}
			else MasterDensityMap[code]+=dens;
		}
	}

	printf("code:\tDensityMap:\tdensityf:\n");
	for (itr=MasterDensityMap.begin();itr!=MasterDensityMap.end();itr++) {
		int code=itr->first;
		double dens=itr->second;
		dens=dens/totvol;

		CresInfo *resinfo=sampler->reslist->GetResInfoPtr(code);
		//double I=2*resinfo->charge-resinfo->baryon-resinfo->strange;
		//real_rhoI+=I*dens;

		printf("%d\t%lf\t%lf\n",code,dens,sampler->densityf[resinfo->ires]);
	}
	//printf("avg_rhoI=%lf\t real_rhoI=%lf\n",avg_rhoI,real_rhoI);
	printf("nparts=%d\n",nparts);
	printf("totvol=%lf\n",totvol);
	printf("nhad=%lf\n",sampler->nhadronsf);
	/****************************************************************************************/

	//Test whether rhoI produced in GetTfMuNH matches rhoI produced in particles (single dummy hyper element)
	/****************************************************************************************
	int nparts=0;
	double avg_rhoI=0.0;
	double real_rhoI=0.0;
	int N=50000*1000;
	list<Chyper *>::iterator it;
	map<int,double>::iterator itr;
	map<int,double> MasterDensityMap;
	double totvol=0.0;

	sampler->totvol=0.0;

	hyper->muB=0.0;
	hyper->T=.080;
	hyper->sigma=93.0;
	hyper->rhoB=0.07;
	hyper->epsilon=(50+1000*hyper->rhoB+750.0*7/15)/1000;
	hyper->rhoI=0.1*hyper->rhoB;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;
	sampler->GetTfMuNH(hyper);
	avg_rhoI=hyper->rhoI;

	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
		totvol+=sampler->totvol;
		for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
			int code=itr->first;
			double dens=itr->second;
			if (MasterDensityMap.count(code)==0) {
				MasterDensityMap.insert(pair<int,double>(code,dens));
			}
			else MasterDensityMap[code]+=dens;
		}
	}

	printf("code:\tDensityMap:\tdensityf0:\n");
	for (itr=MasterDensityMap.begin();itr!=MasterDensityMap.end();itr++) {
		int code=itr->first;
		double dens=itr->second;
		dens=dens/totvol;

		CresInfo *resinfo=sampler->reslist->GetResInfoPtr(code);
		double I=2*resinfo->charge-resinfo->baryon-resinfo->strange;
		real_rhoI+=I*dens;

		printf("%d\t%lf\t%lf\n",code,dens,sampler->densityf0[resinfo->ires]);
	}
	printf("avg_rhoI=%lf\t real_rhoI=%lf\n",avg_rhoI,real_rhoI);
	printf("nparts=%d\n",nparts);
	printf("totvol=%lf\n",sampler->totvol);
	printf("nhad=%lf\n",sampler->nhadronsf0);
	****************************************************************************************/

	//Test viscous corrections
	/*****************************************************************************************
	int nparts=0;
	int N=50000*100;
	map<int,vector<Cpart*>>::iterator it;
	vector<Cpart*>::iterator itr;
	map<int,double>::iterator iter;
	map<int,double> MasterDensityMap;
	double totvol=0.0;


	printf("test 7\n");

	//dummy hyper element
	hyper->tau=12;
	hyper->muB=0.0;
	hyper->T=.137;
	hyper->sigma=93.0;
	hyper->rhoB=0.07;
	hyper->epsilon=.121;
	hyper->rhoI=0.1*hyper->rhoB;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;
	for (int i=0;i<4;i++) hyper->u[i]=0;
	hyper->dOmega[0]=0.089;
	hyper->dOmega[1]=-0.005;
	hyper->dOmega[2]=0.0078;
	hyper->dOmega[3]=0.0;
	hyper->pitilde[0][0]=-0.000167;
	hyper->pitilde[1][1]=-0.000380;
	hyper->pitilde[2][2]=0.000033;
	hyper->pitilde[1][2]=hyper->pitilde[2][1]=0.001156;
	hyper->pitilde[3][3]=-0.000001;
	hyper->pitilde[3][1]=hyper->pitilde[1][3]=0.0;
	hyper->pitilde[3][2]=hyper->pitilde[2][3]=0.0;
	hyper->pitilde[0][1]=hyper->pitilde[1][0]=0.000867;
	hyper->pitilde[0][2]=hyper->pitilde[2][0]=-0.000102;
	hyper->pitilde[0][3]=hyper->pitilde[3][0]=0.0;

	sampler->GetTfMuNH(hyper);
	sampler->CalcLambda();
	printf("test 6\n");

	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
		totvol+=sampler->totvol;
		for (iter=sampler->DensityMap.begin();iter!=sampler->DensityMap.end();iter++) {
			int code=iter->first;
			double dens=iter->second;
			if (MasterDensityMap.count(code)==0) {
				MasterDensityMap.insert(pair<int,double>(code,dens));
			}
			else MasterDensityMap[code]+=dens;
		}
	}
	printf("totvol=%lf\n",totvol);
	printf("test 5\n");
	for (iter=MasterDensityMap.begin();iter!=MasterDensityMap.end();iter++) {
		iter->second=iter->second/totvol;
		printf("code=%d dens=%lf\n",iter->first,iter->second);
	}
	printf("test 4\n");
	double dp=.001;
	double p[4];
	double ptot;
	long double T[4][4];
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			T[i][j]=0.0;
		}
	}
	double E;
	int nsteps=10;
	double avg;
	Cpart *part;
	CresInfo *resinfo;
	double dens;
	int code;
	double temp;
	for (it=sampler->pmap.begin();it!=sampler->pmap.end();it++) {//for each species
		code=it->first;
		if (MasterDensityMap.count(code)==1) {
			dens=MasterDensityMap[code];
		}
		else dens=0.0;
		for (itr=it->second.begin();itr!=it->second.end();itr++) { //for each part of a given species
			part=*itr;
			part->p[0]=(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3])/(2*pow(part->msquared,0.5));
			for (int i=0;i<4;i++) { //row
				for (int j=i;j<4;j++) { //column
					temp=dens*(part->p[i]*part->p[j]/part->p[0])/nparts;
					//printf("dens=%lf p[%d]=%lf p[%d]=%lf p[0]=%lf \n",dens,i,part->p[i],j,part->p[j],part->p[0]);
					T[i][j]+=temp;
				}
			}
		}
	}
	printf("T= \n");
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			printf("%Lf ",T[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("T-P= \n");
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			if (i==j) T[i][j]-=sampler->Pf;
			printf("%Lf ",T[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("nparts=%d\n",nparts);
	*****************************************************************************************/

	//Test GetTfMuNH with real hyper-elements
	/*****************************************************************************************
	list<Chyper *>::iterator itr;
	int n=0;
	int nparts=0;
	for (itr=ms.hyperlist.begin();itr!=ms.hyperlist.end();itr++){
		if (n%500==0) {
			hyper=*itr;
			printf("epsilon=%g, rhoB=%6.4f, nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n",hyper->epsilon,hyper->rhoB,hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
			sampler->GetTfMuNH(hyper);
			nparts+=sampler->MakeParts(hyper);
			printf("epsilon=%g, rhoB=%6.4f, nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n\n",hyper->epsilon,hyper->rhoB,hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
		}
		n++;
	}
	*****************************************************************************************/

	//Test GetTfMuNH with dummy hyper-elements
	/*****************************************************************************************
	for(int ihyper=0;ihyper<=NHYPER;ihyper++){
		hyper->muB=0.0;
		hyper->T=.080;
		hyper->sigma=93.0;
		for(int jhyper=0;jhyper<=NHYPER;jhyper++){
			hyper->rhoB=jhyper*0.15/double(NHYPER);
			hyper->epsilon=(50+1000*hyper->rhoB+750.0*ihyper/double(NHYPER))/1000;
			hyper->rhoI=0.1*hyper->rhoB;
			hyper->rhoS=0.0;
			sampler->GetTfMuNH(hyper);
			sampler->CalcLambda();
			printf("epsilon=%g, rhoB=%6.4f, nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n",hyper->epsilon,hyper->rhoB,hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
			printf("lambdaf=%lf\n",sampler->lambdaf);
		}
		printf("--------------------------------------------------------------------------------------------\n");
	}
	printf("------ hyper elements created ----------- \n");
	//for (int i=0; i<493; i++) {
	//	printf("%d\t%lf\n",i,sampler->densityf0[i]);
	//}
	*****************************************************************************************/

	//This is Scott's, I don't know what it's for.
	/*****************************************************************************************
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	*****************************************************************************************/

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
