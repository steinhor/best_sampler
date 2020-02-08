#include "pratt_sampler/master.h"
#include "pratt_sampler/constants.h"
using namespace std;
using namespace pratt_sampler;
CmeanField *CmasterSampler::meanfield=NULL;

CmasterSampler::CmasterSampler(CparameterMap *parmapin){
	parmap=parmapin;
	randy=new Crandy(-1234);
	reslist=new CresList(parmap);
	NEVENTS=0;
	RESWIDTH_ALPHA=parmap->getD("SAMPLER_RESWIDTH_ALPHA",0.5);
	NTF=parmap->getI("SAMPLER_NTF",1);
	NSIGMAF=parmap->getI("SAMPLER_NSIGMAF",1);
	TFmin=parmap->getD("SAMPLER_TFMIN",155.0);
	TFmax=parmap->getD("SAMPLER_TFMAX",155.0);
	SIGMAFmin=parmap->getD("SAMPLER_SIGMAFMIN",93.0);
	SIGMAFmax=parmap->getD("SAMPLER_SIGMAFMAX",93.0);
	SETMU0=parmap->getB("SAMPLER_SETMU0",false);
	CALCMU=parmap->getB("SAMPLER_CALCMU",false);
	FINDT=parmap->getB("SAMPLER_FINDT", false);
	NEVENTS_TOT=parmap->getI("SAMPLER_NEVENTS_TOT",1);
	NEVENTS=0; // running count of events
	DELTF=(TFmax-TFmin)/double(NTF);
	if(NTF==0)
		DELTF=0.0;
	DELSIGMAF=(SIGMAFmax-SIGMAFmin)/double(NSIGMAF);
	if(NSIGMAF==0)
		DELSIGMAF=0.0;
	MEANFIELD=parmap->getS("SAMPLER_MEANFIELD","simple");
	if(MEANFIELD=="simple")
		meanfield=new CmeanField_Simple(parmap);
	else{
		printf("Don't recognize SAMPLER_MEANFIELD from parameter map=%s\n",MEANFIELD.c_str());
		exit(1);
	}
	Csampler::randy=randy;
	CresInfo::randy=randy;
	Csampler::mastersampler=this;
	Csampler::reslist=reslist;
	Csampler::parmap=parmap;
	Csampler::CALCMU=CALCMU;
	Csampler::bose_corr=parmap->getB("SAMPLER_BOSE_CORR",false);
	Csampler::n_bose_corr=parmap->getI("SAMPLER_N_BOSE_CORR",1);
	int it,isigma;
	hyperlist.clear();
	sampler.resize(NTF+1);
	for(it=0;it<=NTF;it++){
		sampler[it].resize(NSIGMAF+1);
		for(isigma=0;isigma<=NSIGMAF;isigma++){
			//printf("calling Csampler()\n");
			sampler[it][isigma]=new Csampler(TFmin+it*DELTF,SIGMAFmin+isigma*DELSIGMAF);
		}
	}
}

CmasterSampler::~CmasterSampler(){
	for(auto p : hyperlist)
		delete p;
	hyperlist.clear();

	Csampler *sampleri;
	int iT,isigma;
	for(iT=0;iT<NTF;iT++){
		sampler[iT].resize(NSIGMAF);
		for(isigma=0;isigma<NSIGMAF;isigma++){
			//printf("it=%d TFmin=%lf DELTF=%lf TFmin+(it+0.5)*DELTF=%lf\n",it,TFmin,DELTF,TFmin+(it+0.5)*DELTF);
			sampleri=sampler[iT][isigma];
			delete sampleri;
		}
		sampler[iT].clear();
		sampler[iT].shrink_to_fit();
	}
	sampler.clear();
	delete reslist;
	delete randy;
}

int CmasterSampler::MakeEvent(){
	int np,nparts=0;
	Chyper *hyper;
	Csampler *samplerptr;
	double Omega0Sum=0.0;
	partlist->Reset();
	Csampler *samplertemp;
	list<Chyper *>::iterator it;
	if(FINDT && NEVENTS==0){
		samplertemp=new Csampler(150.0,0.093);
	}
	for(it=hyperlist.begin();it!=hyperlist.end();it++){
		hyper=*it;
		Omega0Sum+=hyper->dOmega[0];
		if(NEVENTS==0 && FINDT){
			samplertemp->GetTfMuNH(hyper);
			hyper->T=samplertemp->Tf;
		}
		samplerptr=ChooseSampler(hyper);
		if(samplerptr->FIRSTCALL){
			hyper->T=samplerptr->Tf;
			samplerptr->GetNHMu0();
			samplerptr->CalcDensitiesMu0();
			samplerptr->CalcLambdaMu0();
			samplerptr->FIRSTCALL=false;
		}
		if(NEVENTS==0){
			if(!SETMU0 && CALCMU){
				samplerptr->GetMuNH(hyper);
				samplerptr->GetNHadronsEpsilonPF(hyper);
				hyper->lambda=samplerptr->CalcLambdaF(hyper);
			}
			if(!SETMU0 && !CALCMU){
				samplerptr->GetNHadronsEpsilonPF(hyper);
				hyper->lambda=samplerptr->CalcLambdaF(hyper);
			}
			if(SETMU0){
				hyper->lambda=samplerptr->lambda0;
			}
		}
		np=samplerptr->MakeParts(hyper);
		nparts+=np;
	}
	NEVENTS+=1;
	if(FINDT && NEVENTS==0){
		delete samplertemp;
	}
	return nparts;
}

Csampler* CmasterSampler::ChooseSampler(Chyper *hyper){
	double T,sigma,del;
	int it,isigma;
	if(NTF==0){
		it=0;
	}
	else{
		T=hyper->T;
		it=floorl((T-TFmin)/DELTF);
		del=T-(TFmin+it*DELTF);
		if(randy->ran()<del/DELTF)
			it+=1;
		if(it<0)
			it=0;
		if(it>NTF)
			it=NTF;
	}
	if(NSIGMAF==0){
		isigma=0;
	}
	else{
		sigma=hyper->sigma;
		isigma=floorl((sigma-SIGMAFmin)/DELSIGMAF);
		del=sigma-(SIGMAFmin+isigma*DELSIGMAF);
		if(randy->ran()<del/DELSIGMAF)
			isigma+=1;
		if(isigma<0)
			isigma=0;
		if(isigma>=NSIGMAF)
			isigma=NSIGMAF-1;
	}
	return sampler[it][isigma];
}

void CmasterSampler::MakeDummyHyper(){
	Chyper *hyper=new Chyper();
	hyperlist.push_back(hyper);
}

void CmasterSampler::ReadHyper2D(){
	string filename;
	Chyper *elem;
	double PIbulk;
	int ielement=0;
	double u0,ux,uy,x,y,udotdOmega,dOmega0,dOmegaX,dOmegaY,pitildexx,pitildeyy,pitildexy,tau,epsilonf;
	double eta,dOmegaZ,uz,Tdec,muB,muS,muC,Pdec;
	double pitilde00,pitilde0x,pitilde0y,pitilde0z,pitildexz,pitildeyz,pitildezz;
	double qmu0,qmu1,qmu2,qmu3;
	double rhoB;

	nelements=0;
	filename=parmap->getS("HYPER_INFO_FILE",string("../local/include/surface_2D.dat"));
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"rb");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open hyper info file\n");
		printf("Error %d \n", errno);
		exit(1);
	}

	double TotalVolume=0.0;

	while(!feof(fptr)){
		elem=new Chyper();
		//read from binary file
		float array[34];
		fread(array, sizeof(array),1,fptr);

		tau = array[0];
		x = array[1];
		y = array[2];
		eta = array[3];
		dOmega0 = array[4];
		dOmegaX = array[5];
		dOmegaY = array[6];
		dOmegaZ = array[7];
		u0 = array[8];
		ux = array[9];
		uy = array[10];
		uz = array[11]/tau;
		u0=sqrt(1.0+ux*ux+uy*uy+uz*uz);

		epsilonf = array[12]*HBARC; //was labeled Edec--guessed this was epsilon
		Tdec = array[13]*HBARC;
		muB  = array[14]*HBARC;
		muS  = array[15]*HBARC;
		muC  = array[16]*HBARC;
		Pdec = array[17]*Tdec - epsilonf;

		pitilde00 = array[18]*HBARC;  // GeV/fm^3
		pitilde0x = array[19]*HBARC;  // GeV/fm^3
		pitilde0y = array[20]*HBARC;  // GeV/fm^3
		pitilde0z = array[21]*HBARC/tau;  // GeV/fm^3
		pitildexx = array[22]*HBARC;  // GeV/fm^3
		pitildexy = array[23]*HBARC;  // GeV/fm^3
		pitildexz = array[24]*HBARC/tau;  // GeV/fm^3
		pitildeyy = array[25]*HBARC;  // GeV/fm^3
		pitildeyz = array[26]*HBARC/tau;  // GeV/fm^3
		pitildezz = array[27]*HBARC/(tau*tau);  // GeV/fm^3

		PIbulk = array[28]*HBARC;   // GeV/fm^3
		rhoB = array[29];  // 1/fm^3

		qmu0 = array[30];
		qmu1 = array[31];
		qmu2 = array[32];
		qmu3 = array[33];

		udotdOmega=tau*(u0*dOmega0+ux*dOmegaX+uy*dOmegaY+uz*dOmegaZ);

		if(!(udotdOmega<0.0)) {
			elem->tau=tau;
			elem->dOmega[0]=dOmega0; //*2.0*b3d->ETAMAX;
			elem->dOmega[1]=dOmegaX; //*2.0*b3d->ETAMAX;
			elem->dOmega[2]=dOmegaY; //*2.0*b3d->ETAMAX;
			elem->dOmega[3]=dOmegaZ; //*2.0*b3d->ETAMAX;

			elem->udotdOmega=udotdOmega; //*2.0*b3d->ETAMAX;

			elem->r[1]=x;
			elem->r[2]=y;
			elem->u[0]=sqrt(1.0+ux*ux+uy*uy+uz*uz);
			elem->u[1]=ux;
			elem->u[2]=uy;
			elem->u[3]=uz;

			elem->muB=muB;
			elem->muS=muS;

			elem->pitilde[0][0]=pitilde00;
			elem->pitilde[1][1]=pitildexx;
			elem->pitilde[2][2]=pitildeyy;
			elem->pitilde[1][2]=elem->pitilde[2][1]=pitildexy;
			elem->pitilde[3][3]=-pitildezz;
			elem->pitilde[3][1]=elem->pitilde[1][3]=pitildexz;
			elem->pitilde[3][2]=elem->pitilde[2][3]=pitildeyz;
			elem->pitilde[0][1]=elem->pitilde[1][0]=pitilde0x;
			elem->pitilde[0][2]=elem->pitilde[2][0]=pitilde0y;
			elem->pitilde[0][3]=elem->pitilde[3][0]=pitilde0z;
			elem->epsilon=epsilonf;
			elem->T=Tdec;
			elem->rhoB=rhoB;
			elem->rhoS=0.0;

			elem->qmu[0]=qmu0;
			elem->qmu[1]=qmu1;
			elem->qmu[2]=qmu2;
			elem->qmu[3]=qmu3;
			elem->rhoB=rhoB;

			hyperlist.push_back(elem);
			TotalVolume+=udotdOmega;
			ielement+=1;

		}
	}
	nelements=ielement;
	printf("Exiting ReadHyper2D() happily, TotalVolume=%lf, nelements=%d\n",TotalVolume,nelements);
}
