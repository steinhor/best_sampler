#include "test.h"
#include "sampler.h"

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
  printf("nhad=%lf\n",sampler->nhadronsf);

}

void test::n_dummy(CmasterSampler ms) {
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

	printf("code: dens: densityf: ires:\n");
	for (itr=sampler->DensityMap.begin();itr!=sampler->DensityMap.end();itr++) {
		itr->second=itr->second/sampler->totvol;
		resinfo=sampler->reslist->GetResInfoPtr(itr->first);
		printf("%d %lf %lf %d\n",itr->first,itr->second,sampler->densityf[resinfo->ires],resinfo->ires);
	}

	printf("nparts=%d\t",nparts);
	printf("totvol=%lf\t",sampler->totvol);
	printf("nparts/totvol=%lf\t",nparts/sampler->totvol);
	printf("nhad=%lf\n",sampler->nhadronsf);

}

void test::bose_terms() {
  Csampler *sampler=new Csampler(.150,0.093);
  Chyper *hyper=new Chyper();

  int N=50000*1000;
	int nparts=0;
	double nhadronsf0,totvol=0.0;
	map<int,double>::iterator iter;
	map<int,double> MasterDensityMap;
	map<double,double>::iterator itr;
  vector<double> dN_dp_p2;

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
			dN_dp_p2.push_back(0.0);
		}

		for (int i=0;i<N;i++) {
			nparts+=sampler->MakeParts(hyper);
      
		}
		//sampler->totvol=N/.180011;
		//for (itr=sampler->dpmap.begin();itr!=sampler->dpmap.end();itr++) {
		//	printf("%lf %lf\n",itr->first,itr->second/sampler->totvol);
		//}

		for (int ibin=0;ibin<ndp;ibin++) {
			printf("%lf %lf\n",ibin*dp+dp/2,dN_dp_p2[ibin]/sampler->totvol);
		}

		printf("nparts=%d totvol=%lf dens=%lf nhadronsf=%lf\n",nparts,sampler->totvol,nparts/sampler->totvol,sampler->nhadronsf);
		//printf("%lf %lf\n",hyper->T,nhadronsf0);
	//}
}
