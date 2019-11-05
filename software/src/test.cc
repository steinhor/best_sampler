#include "pratt_sampler/test.h"
#include "pratt_sampler/sampler.h"

void test::n_hyperlist(CmasterSampler ms,int nsamplers) {
  Csampler *sampler=new Csampler(0.150,0.0930);
  Chyper *hyper=*ms.hyperlist.begin();

  int nparts=0;
  int temp;
  int n=0;
  list<Chyper *>::iterator it;
  map<int,double>::iterator itr;
  map<int,double> MasterDensityMap;
  double totvol=0.0;

  for (int i=0;i<nsamplers;i++) {
    if (i!=0) {
      sampler=new Csampler(0.135,0.0930);
      sampler->partmap=new CpartMap;
    }
    sampler->totvol=0.0;
    sampler->GetTfMuNH(hyper);
    for (it=ms.hyperlist.begin();it!=ms.hyperlist.end();it++) { //for each hyper element
      hyper=*it;
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
      else {
        MasterDensityMap[code]+=dens;
      }
    }
  }

  printf("code:\tDensityMap:\tdensityf:\n");
  for (itr=MasterDensityMap.begin();itr!=MasterDensityMap.end();itr++) {
    int code=itr->first;
    double dens=itr->second;
    dens=dens/totvol;

    CresInfo *resinfo=sampler->reslist->GetResInfoPtr(code);
    printf("%d\t%lf\t%lf\n",code,dens,sampler->densityf[resinfo->ires]);
  }
  printf("nparts=%d\n",nparts);
  printf("totvol=%lf\n",totvol);
	printf("nparts/totvol=%lf\t",nparts/totvol);
  printf("nhadronsf=%lf\n",sampler->nhadronsf);
  printf("nhadronsf0=%lf\n",sampler->nhadronsf0);

}

void test::n_dummy(bool readout) {
  Csampler *sampler=new Csampler(0.150,0.0930);
  Chyper *hyper=new Chyper();

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

  if(readout) {
  	printf("code: dens: densityf: ires:\n");
  	for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
  		itr->second=itr->second/sampler->totvol;
  		resinfo=sampler->reslist->GetResInfoPtr(itr->first);
  		printf("%d %lf %lf %d\n",itr->first,itr->second,sampler->densityf[resinfo->ires],resinfo->ires);
  	}
  }

	printf("nparts=%d\t",nparts);
	printf("totvol=%lf\t",sampler->totvol);
	printf("nparts/totvol=%lf\t",nparts/sampler->totvol);
  printf("nhadronsf=%lf\n",sampler->nhadronsf);
  printf("nhadronsf0=%lf\n",sampler->nhadronsf0);

}

void test::bose_terms() {
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

  FILE *file=fopen("corrections.txt","w");
  calc_bose(sampler,hyper,file);
  fclose(file);

  sampler->bose_test_off=true;

  file=fopen("no_corrections.txt","w");
  calc_bose(sampler,hyper,file);
  fclose(file);
}

void test::calc_bose(Csampler *sampler,Chyper *hyper, FILE *file) {
  int N=50000*1000;
  int nparts=0;
  double nhadronsf0,totvol=0.0;
  vector<double> dN_dp_p2;

	sampler->GetTfMuNH(hyper);
  bool bose_off=sampler->bose_test_off;
	sampler=new Csampler(hyper->T,hyper->sigma);
  sampler->bose_test=true;
  if (bose_off) sampler->bose_test_off=true;
	sampler->GetTfMuNH(hyper);
	sampler->CalcLambdaF();

	printf("T=%lf\n",sampler->Tf);
	nhadronsf0=sampler->nhadronsf0;
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

	printf("nparts=%d totvol=%lf dens=%lf nhadronsf=%lf nhadronsf0=%lf\n",nparts,sampler->totvol,nparts/sampler->totvol,sampler->nhadronsf,sampler->nhadronsf0);

}

void test::test_viscous() {
  int nparts=0;
	int N=50000*10;
	map<int,vector<Cpart*>>::iterator it;
	vector<Cpart*>::iterator itr;
	double totvol=0.0;

  double p[4];
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
	hyper->pitilde[1][2]=hyper->pitilde[2][1]=.005;
	hyper->pitilde[3][3]=0.0;
	hyper->pitilde[3][1]=hyper->pitilde[1][3]=0.0;
	hyper->pitilde[3][2]=hyper->pitilde[2][3]=0.0;
	hyper->pitilde[0][1]=hyper->pitilde[1][0]=0.0;
	hyper->pitilde[0][2]=hyper->pitilde[2][0]=0.0;
	hyper->pitilde[0][3]=hyper->pitilde[3][0]=0.0;

	Csampler *sampler=new Csampler(hyper->T,hyper->sigma);
  sampler->viscous_test=true;

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

}

void test::newton_method_hyperlist(CmasterSampler ms) {
  Csampler *sampler=new Csampler(0.150,0.0930);
  Chyper *hyper=new Chyper();
  list<Chyper *>::iterator itr;
	int n=0;
	int nparts=0;
	for (itr=ms.hyperlist.begin();itr!=ms.hyperlist.end();itr++){
		if (n%500==0) {
			hyper=*itr;
			printf("epsilon=%g, rhoB=%6.4f, rhoI=%6.4f, rhoS=%6.4f, T=%7.3f\n",hyper->epsilon,hyper->rhoB,hyper->rhoI,hyper->rhoS,hyper->T);
			sampler->GetTfMuNH(hyper);
			printf("nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n\n",hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
		}
		n++;
	}
}

void test::newton_method_dummy() {
  Csampler *sampler=new Csampler(0.150,0.0930);
  Chyper *hyper=new Chyper();
  const int NHYPER=15;

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
}

void test::eps_rho_derivs() {
  Csampler *sampler=new Csampler(0.150,0.0930);

  Eigen::MatrixXd A(4,4);
	Eigen::MatrixXd B(4,4);
	Eigen::VectorXd dmu(4);
	dmu[0]=0.05; dmu[1]=dmu[2]=dmu[3]=0.002;
	Eigen::VectorXd mu(4);
	Eigen::VectorXd rho(4),rho1(4),rho2(4);
	double rhoB,rhoI,rhoS,epsilon;

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
}
