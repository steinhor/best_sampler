#include "pratt_sampler/master.h"
using namespace std;

int main(int argc, char *argv[]){
	int nprint=100,iprint=0;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	Csampler *sampler=new Csampler(0.150,0.0930);
	sampler->partmap=new CpartMap;
	sampler->randy=new Crandy(1234);
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
	/****************************************************************************************
	int nparts=0;
	int temp;
	int n=0;
	double avg_rhoI=0.0;
	double real_rhoI=0.0;
	list<Chyper *>::iterator it;
	map<int,double>::iterator itr;
	map<int,double> MasterDensityMap;
	double totvol=0.0;

	for (int i=0;i<1;i++) { //for each sampler
		if (i!=0) {
			sampler=new Csampler(0.135,0.0930);
			sampler->partmap=new CpartMap;
		}
		sampler->totvol=0.0;
		for (it=ms.hyperlist.begin();it!=ms.hyperlist.end();it++) { //for each hyper element
			hyper=*it;
			sampler->GetTfMuNH(hyper);
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
	****************************************************************************************/

	//Test densities with fake hyper element
	/****************************************************************************************/
	int nparts=0;
	int N=50000*100;
	map<int,double>::iterator itr;
	CresInfo *resinfo;

	sampler->totvol=0.0;
	hyper->muB=0.0;
	hyper->T=.140;
	hyper->sigma=93.0;
	hyper->rhoB=0.07;
	hyper->epsilon=(50+1000*hyper->rhoB+750.0*7/15)/1000;
	hyper->rhoI=0.1*hyper->rhoB;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;
	sampler->GetTfMuNH(hyper);

	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
	}

	//printf("code: dens: densityf: ires:\n");
	//for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
	//	itr->second=itr->second/sampler->totvol;
	//	resinfo=sampler->reslist->GetResInfoPtr(itr->first);
	//	printf("%d %lf %lf %d\n",itr->first,itr->second,sampler->densityf[resinfo->ires],resinfo->ires);
	//}

	printf("nparts=%d\t",nparts);
	printf("totvol=%lf\t",sampler->totvol);
	printf("nparts/totvol=%lf\t",nparts/sampler->totvol);
	printf("nhad=%lf\n",sampler->nhadronsf);
	/****************************************************************************************/

	//produce data for graph to test bose corrections
	/*****************************************************************************************
	int N=50000*1000;
	int nparts=0;
	double nhadronsf0,totvol=0.0;
	map<int,double>::iterator iter;
	map<int,double> MasterDensityMap;
	map<double,double>::iterator itr;

	hyper->muB=0.0;
	hyper->sigma=.0930;
	hyper->rhoB=0.007;
	hyper->epsilon=.121;
	hyper->rhoI=0.1*hyper->rhoB;
	hyper->rhoS=0.0;
	hyper->udotdOmega=0.2;

	//for(int ihyper=0;ihyper<=NHYPER;ihyper++){ //make NHYPER elements, each w/ a diff temp
		//printf("starting element %d\t",ihyper);
		//hyper->T=.1*(1+ihyper/double(NHYPER));
		hyper->T=.150;
		sampler->GetTfMuNH(hyper);
		sampler=new Csampler(hyper->T,hyper->sigma);
		sampler->GetTfMuNH(hyper);
		sampler->CalcLambdaF();

		printf("T=%lf\n",sampler->Tf);
		nhadronsf0=sampler->nhadronsf0;
		MasterDensityMap.clear();
		nparts=0;

		double ndp=400,pbound=4.5;
		double dp=pbound/ndp,p,Tf=.144535,m=.138;
		int ibin;
		Cpart *part;
		for(int i=0;i<ndp;i++) {
			sampler->dN_dp_p2.push_back(0.0);
		}

		for (int i=0;i<N;i++) {
			nparts+=sampler->MakeParts(hyper);
		}
		//sampler->totvol=N/.180011;
		//for (itr=sampler->dpmap.begin();itr!=sampler->dpmap.end();itr++) {
		//	printf("%lf %lf\n",itr->first,itr->second/sampler->totvol);
		//}

		for (int ibin=0;ibin<ndp;ibin++) {
			printf("%lf %lf\n",ibin*dp+dp/2,sampler->dN_dp_p2[ibin]/sampler->totvol);
		}

		printf("nparts=%d totvol=%lf dens=%lf nhadronsf=%lf\n",nparts,sampler->totvol,nparts/sampler->totvol,sampler->nhadronsf);
		//printf("%lf %lf\n",hyper->T,nhadronsf0);
	//}
	*****************************************************************************************/

	//Test viscous corrections
	/*****************************************************************************************
	int nparts=0;
	int N=50000*10;
	map<int,vector<Cpart*>>::iterator it;
	vector<Cpart*>::iterator itr;
	map<int,double>::iterator iter;
	map<double,double>::iterator itt;
	map<int,double> MasterDensityMap;
	double totvol=0.0;
	CresInfo *resinfo;

	//dummy hyper element
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
	hyper->pitilde[1][2]=hyper->pitilde[2][1]=0.0005;
	hyper->pitilde[3][3]=0.0;
	hyper->pitilde[3][1]=hyper->pitilde[1][3]=0.0;
	hyper->pitilde[3][2]=hyper->pitilde[2][3]=0.0;
	hyper->pitilde[0][1]=hyper->pitilde[1][0]=0.0;
	hyper->pitilde[0][2]=hyper->pitilde[2][0]=0.0;
	hyper->pitilde[0][3]=hyper->pitilde[3][0]=0.0;
	sampler=new Csampler(hyper->T,hyper->sigma);
	for (int j=0;j<1;j++) {
		hyper->rhoB=.007;
		hyper->rhoI=0.1*hyper->rhoB;
		//hyper->rhoB=.000007*pow(10.0,4.0*double(j)/100.0);
		sampler->GetTfMuNH(hyper);
		printf("hyper->T=%lf\tsampler->Tf=%lf\n",hyper->T,sampler->Tf);
		sampler->CalcLambda();
		//printf("%lf ",hyper->rhoB);
		for (int i=0;i<N;i++) {
			hyper->ihyp=i;
			nparts+=sampler->MakeParts(hyper);
		}
	}
	totvol=sampler->totvol;
	for (itt=sampler->dpmap.begin();itt!=sampler->dpmap.end();itt++) {
		printf("%lf %lf\n",itt->first,itt->second*(2*PI*PI*pow(HBARC,3))/totvol);
	}
	for (iter=sampler->DensityMap.begin();iter!=sampler->DensityMap.end();iter++) {
		iter->second=iter->second/totvol;
		//resinfo=sampler->reslist->GetResInfoPtr(iter->first);
		//printf("code=%d dens=%lf densityf=%lf ires=%d\n",iter->first,iter->second,sampler->densityf[resinfo->ires],resinfo->ires);
	}
	double p[4];
	double ptot;
	double T[4][4];
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			T[i][j]=0.0;
		}
	}
	Cpart *part;
	double dens;
	int code;
	double temp;
	double psquared;
	double P;
	double Tcalc=0;
	bool x=false;


	for (it=sampler->pmap.begin();it!=sampler->pmap.end();it++) {//for each species
		code=it->first;
		if (sampler->DensityMap.count(code)==1) {
			dens=sampler->DensityMap[code];
		}
		else dens=0.0;
		for (itr=it->second.begin();itr!=it->second.end();itr++) { //for each part of a given species
			part=*itr;
			part->p[0]=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3]+part->msquared);
			psquared=part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3];
			if (!isnan(part->p[0])) P+=sampler->nhadronsf*(psquared/(3*part->p[0]))/(nparts);
			if (!isnan(part->p[0])) Tcalc+=(psquared/(3*part->p[0]))/(nparts);
			else {
				printf("psquared=%lf part->p[0]=%lf nparts=%d\n",psquared,part->p[0],nparts);
			}
			for (int i=0;i<4;i++) { //row
				for (int j=i;j<4;j++) { //column
					temp=sampler->nhadronsf*(part->p[i]*part->p[j]/part->p[0])/nparts;
					if (!isnan(part->p[0])) T[i][j]+=temp;
				}
			}
		}
	}
	double pbound=2.5;
	double ndp=100.0;
	resinfo=sampler->reslist->GetResInfoPtr(321);
	double m=resinfo->mass;
	double Tf=0.144;
	int alpha,beta;
	FourVector pnovisc;
	for (long int i=0;i<N;i++) {
		part=new Cpart;
		sampler->randy->generate_boltzmann(m,Tf,part->p);
		for (double i=0;i<=pbound;i+=pbound/double(ndp)) {
			if (part->p[0]>=i && part->p[0]<(i+pbound/ndp)) {
				if (sampler->dpmap.count(i+pbound/(2*ndp))==0) sampler->dpmap.insert(pair<double,double>(i+pbound/(2*ndp),1));
				else sampler->dpmap[i+pbound/(2*ndp)]+=1;
				//printf("binned! i=%lf\n",i);
				break;
			}
		}
        //part->p[0]=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3]+m*m);
		//printf("p0calc=%lf\n",sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3]+m*m));

		psquared=part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3];
		//P+=sampler->nhadronsf*(psquared/(3*part->p[0]))/N;
		Tcalc+=(psquared/(3*part->p[0]))/N;
		for (int i=0;i<4;i++) { //row
			for (int j=i;j<4;j++) { //column
				temp=(part->p[i]*part->p[j]/part->p[0])/N; //sampler->nhadronsf
				//printf("dens=%lf p[%d]=%lf p[%d]=%lf p[0]=%lf \n",dens,i,part->p[i],j,part->p[j],part->p[0]);
				T[i][j]+=temp;
			}
		}
	}
	for (map<double,double>::iterator it=sampler->dpmap.begin();it!=sampler->dpmap.end();it++) {
		//printf("%lf %lf\n",it->first,it->second);
	}
	P=Tcalc;

	printf("T= \n");
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			printf("%lf ",T[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("P=%lf\tsampler->Pf=%lf\tTcalc=%lf\tTsampler=%lf\n\n",P,sampler->Pf,Tcalc,sampler->Pf/sampler->nhadronsf);
	printf("\n");
	printf("T-Pcalc= \n");
	for (int i=0;i<4;i++) {
		for (int j=0;j<4;j++) {
			if (i==j) T[i][j]-=P;
			printf("%lf ",T[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("totvol=%lf densityf=%lf\n",totvol,sampler->nhadronsf);
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
