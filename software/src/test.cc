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

  printf("code:\tDensityMap:\tdensity0:\n");
  for (itr=MasterDensityMap.begin();itr!=MasterDensityMap.end();itr++) {
    int code=itr->first;
    double dens=itr->second;
    dens=dens/totvol;

    CresInfo *resinfo=sampler->reslist->GetResInfoPtr(code);
    printf("%d\t%lf\t%lf\n",code,dens,sampler->density0i[resinfo->ires]);
  }
  printf("nparts=%d\n",nparts);
  printf("totvol=%lf\n",totvol);
	printf("nparts/totvol=%lf\t",nparts/totvol);
  //printf("nhadronsf=%lf\n",sampler->nhadrons);
  printf("nhadrons0=%lf\n",sampler->nhadrons0);

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
  sampler->GetNHadronsEpsilonPF(hyper);
  sampler->CalcLambdaF(hyper);
  sampler->CalcDensitiesMu0();

	for (int i=0;i<N;i++) {
		nparts+=sampler->MakeParts(hyper);
	}

  if(readout) {
  	printf("code: dens: density0i: ires:\n");
  	for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
  		itr->second=itr->second/sampler->totvol;
  		resinfo=sampler->reslist->GetResInfoPtr(itr->first);
  		printf("%d %lf %lf %d\n",itr->first,itr->second,sampler->density0i[resinfo->ires],resinfo->ires);
  	}
  }

	printf("nparts=%d\t",nparts);
	printf("totvol=%lf\t",sampler->totvol);
	printf("nparts/totvol=%lf\t",nparts/sampler->totvol);
  //printf("nhadrons=%lf\n",sampler->nhadrons);
  printf("nhadrons0=%lf\n",sampler->nhadrons0);

}

void test::bose_terms(CmasterSampler *ms) {
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
  calc_bose(ms,sampler,hyper,file);
  fclose(file);

  sampler->bose_test_off=true;

  file=fopen("no_corrections.txt","w");
  calc_bose(ms,sampler,hyper,file);
  fclose(file);
}

void test::calc_bose(CmasterSampler *ms,Csampler *sampler,Chyper *hyper, FILE *file) {
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
			sampler->CalcLambdaF(hyper);
			printf("epsilon=%g, rhoB=%6.4f, nhadrons=%7.4f, T=%7.3f, muB=%7.4f, muI=%7.4f, muS=%7.4f \n",hyper->epsilon,hyper->rhoB,hyper->nhadrons,hyper->T,hyper->muB,hyper->muI,hyper->muS);
			printf("lambda0=%lf\n",sampler->lambda0);
		}
		printf("--------------------------------------------------------------------------------------------\n");
	}
	printf("------ hyper elements created ----------- \n");
	//for (int i=0; i<493; i++) {
	//	printf("%d\t%lf\n",i,sampler->densityf0[i]);
	//}
}
