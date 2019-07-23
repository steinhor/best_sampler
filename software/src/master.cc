#include "master.h"
using namespace std;
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
	if(NTF==1)
		DELTF=0.0;
	else
		DELTF=(TFmax-TFmin)/double(NTF);
	if(NSIGMAF==1)
		DELSIGMAF=0.0;
	else
		DELSIGMAF=(SIGMAFmax-SIGMAFmin)/double(NSIGMAF);
	Csampler::randy=randy;
	Csampler::mastersampler=this;
	Csampler::reslist=reslist;
	Csampler::parmap=parmap;
	Cpart::reslist=reslist;
	int it,isigma;
	sampler.resize(NTF);
	for(it=0;it<NTF;it++){
		sampler[it].resize(NSIGMAF);
		for(isigma=0;isigma<NSIGMAF;isigma++){
			sampler[it][isigma]=new Csampler(TFmin+(it+0.5)*DELTF,SIGMAFmin+isigma*DELSIGMAF);
		}
	}
}

CmasterSampler::~CmasterSampler(){
	//printf("MasterSampler ending, NEVENTS=%d\n",NEVENTS);
}

int CmasterSampler::MakeEvent(){
	int np;
	Csampler *sampler;
	double Omega0Sum=0.0;
	int nelements=hyperlist.size();
	int ielement,nparts=0;
	list<Chyper *>::iterator it;
	for(it=hyperlist.begin();it!=hyperlist.end();it++){
		Omega0Sum+=(*it)->dOmega[0];
		np=0;
		sampler=ChooseSampler(*it);
		if(!SETMU0){
			if(NEVENTS==0){
				sampler->GetTfMuNH(*it);
				sampler->CalcLambdaF();
				(*it)->muB=sampler->muB;
				(*it)->muI=sampler->muI;
				(*it)->muS=sampler->muS;
				(*it)->nhadrons=sampler->nhadronsf;
				(*it)->epsilon=sampler->epsilonf;
				(*it)->lambda=sampler->lambdaf;
			}
			else{
				sampler->muB=(*it)->muB;
				sampler->muI=(*it)->muI;
				sampler->muS=(*it)->muS;
				sampler->nhadronsf=(*it)->nhadrons;
				sampler->epsilonf=(*it)->epsilon;
				sampler->lambdaf=(*it)->lambda;
				sampler->Pf=sampler->nhadronsf*sampler->Tf;
			}
		}
		np+=sampler->MakeParts(*it);
		nparts+=np;
	}
	NEVENTS+=1;
	return nparts;
}

Csampler* CmasterSampler::ChooseSampler(Chyper *hyper){
	double T,sigma,rhoB,rhoI,T0,T1;
	int it,isigma;
	if(NTF==1){
		it=0;
	}
	else{
		T=hyper->T;
		it=lrint((T-TFmin-0.5*DELTF)/DELTF);
		if(it<0)
			it=0;
		if(it>=NTF)
			it=NTF-1;
	}
	if(NSIGMAF==1){
		isigma=0;
	}
	else{
		sigma=hyper->sigma;
		isigma=lrint((sigma-SIGMAFmin-0.5*DELSIGMAF)/DELSIGMAF);
		if(isigma<0)
			isigma=0;
		if(isigma>=NSIGMAF)
			sigma=NSIGMAF-1;
	}
	return sampler[it][isigma];
}

/*
void CmasterSampler::CalcPiFromParts(){
	double pi[4][4]={0.0};
	CpartMap::iterator ppos;
	CPart *part;
	double *p,pressure,Ptest=0.0,etest=0.0;
	int ipart,alpha,beta;
	double volume=nelements*hyper[0].Omega[0];
	for(ppos=b3d->DeadPartMap.begin();ppos!=b3d->PartMap.end();ppos++){
		part=ppos->second;
		pressure=(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3])/(3.0*part->p[0]);
		Ptest+=pressure;
		etest+=p[0];
		for(alpha=1;alpha<4;alpha++){
			pi[alpha][alpha]-=pressure;
			for(beta=1;beta<4;beta++){
				pi[alpha][beta]+=p[alpha]*p[beta]/p[0];
			}
		}
	}
	printf("----From particle momenta  -----\n");
	printf("epsilon=%10.4f, P=%10.4f\n",etest/volume,Ptest/volume);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pi[alpha][beta]=pi[alpha][beta]/volume;
			printf("%9.6f ",pi[alpha][beta]/hyper[0].P);
		}
		printf("\n");
	}
	printf("---- From pi_ij  input -------\n");
	printf("epsilon=%10.4f, P=%10.4f\n",hyper[0].epsilon,hyper[0].P);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%9.6f ",hyper[0].pitilde[alpha][beta]/hyper[0].P);
		}
		printf("\n");
	}
}
*/

void CmasterSampler::ReadHyper2D(){
	string filename;
	Chyper *elem;
	double dumbo,udotn,PIbulk;
	double u0,ux,uy,x,y,udotdOmega,dOmega0,dOmegaX,dOmegaY,dOmegaMax,pitildexx,pitildeyy,pitildexy,tau,epsilonf;
	int alpha,beta;
	double sigma,PI;
	int ielement,initarraysize=1000;
	char dummy[300];
	double eta,dOmegaZ,uz,Edec,Tdec,muB,muS,muC,Pdec;
	double pitilde00,pitilde0x,pitilde0y,pitilde0z,pitildexz,pitildeyz,pitildezz;
	double qmu0,qmu1,qmu2,qmu3;
	double rhoB;
	double nproton=parmap->getD("N_PROTON",0.4);

	//element.clear();
	nelements=0;
	//b3d->TotalVolume=0.0;
	filename=parmap->getS("HYPER_INFO_FILE",string("../local/include/surface_2D.dat"));
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"rb");

	ielement=0;
	double TotalVolume=0.0;

	int n=0;
	while(!feof(fptr)){
		elem=new Chyper();
		elem->ihyp=ielement;
		//if(element.size()==ielement) element.resize(element.size()+initarraysize);

		//read from binary file
		float array[34];
		fread(array, sizeof(array),1,fptr);
		/*
		for (int i = 0; i < 34; i++) {
			float temp = 0.;
			fptr.read((char*)&temp, sizeof(float));
			array[i] = temp;
		}
		*/

		//elem=&element[ielement];
		//elem->T=b3d->parmap.getD("FREEZEOUT_TEMP",0.155);

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

		/* OLD LINES FOR READING DATA
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tau,&x,&y,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,&udotdOmega,&dOmegaMax,
		&pitildexx,&pitildeyy,&pitildexy);
		u0=sqrt(1.0+ux*ux+uy*uy);
		*/

		//printf("udotdOmega=%g =? %g\n",udotdOmega,u0*dOmega0+ux*dOmegaX+uy*dOmegaY);

		if(udotdOmega<0.0){
			/*
			printf("udotdOmega<0!!!\n");
			printf("ielement=%d, u=(%lf,%lf,%lf,%lf), dOmega=(%lf,%lf,%lf,%lf)\n",ielement, u0, ux, uy, uz, dOmega0, dOmegaX, dOmegaY, dOmegaZ);
			//exit(1);
			*/
			n++;
		}
		else {
			elem->tau=tau;
			elem->dOmega[0]=dOmega0; //*2.0*b3d->ETAMAX;
			elem->dOmega[1]=dOmegaX; //*2.0*b3d->ETAMAX;
			elem->dOmega[2]=dOmegaY; //*2.0*b3d->ETAMAX;
			elem->dOmega[3]=dOmegaZ; //*2.0*b3d->ETAMAX;

			elem->udotdOmega=udotdOmega; //*2.0*b3d->ETAMAX;

			elem->r[1]=x;
			elem->r[2]=y;
			elem->u[0]=u0;
			elem->u[1]=ux;
			elem->u[2]=uy;
			elem->u[3]=uz;

			elem->muB=muB;
			elem->muS=muS;

			//elem->CalcOmegaMax();
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
			//elem->density=&densityf;
			elem->T=Tdec; //was Tf
			//elem->P=Pdec; //was Pf
			//elem->lambda=lambdaf;
			elem->rhoB=rhoB;
			elem->rhoS=0.0;
			//elem->rhoI=0.5*(nproton-(1-nproton));

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
	printf("Exiting ReadHyper2D() happily, TotalVolume=%lf, nelements=%d, n=%d\n",TotalVolume,nelements,n);
}
