#include "msu_sampler/master.h"
#include "msu_sampler/constants.h"

//#define __TEST_PITILDE_TRACE__

using namespace std;
using namespace msu_sampler;
CmeanField *CmasterSampler::meanfield=NULL;

CmasterSampler::CmasterSampler(CparameterMap *parmapin){
	parmap=parmapin;
	randy=new Crandy(-1234);
	reslist=new CresList(parmap);
	NEVENTS=0;
	RESWIDTH_ALPHA=parmap->getD("SAMPLER_RESWIDTH_ALPHA",0.5);
	TFmax=0.250;
	DELTF=parmap->getD("SAMPLER_DELTF",0.001);
	NTF=lrint(TFmax/DELTF);
	NSIGMAF=parmap->getD("SAMPLER_NSIGMAF",1);
	if(NSIGMAF<1)
		NSIGMAF=1;
	SIGMAFmin=parmap->getD("SAMPLER_SIGMAFMIN",93.0);
	SIGMAFmax=parmap->getD("SAMPLER_SIGMAFMAX",93.0);
	SETMU0=parmap->getB("SAMPLER_SETMU0",false);
	CALCMU=parmap->getB("SAMPLER_CALCMU",false);
	FINDT=parmap->getB("SAMPLER_FINDT", false);
	NEVENTS_TOT=parmap->getI("SAMPLER_NEVENTS_TOT",1);
	NEVENTS=0; // running count of events
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
	Csampler::BJORKEN_2D=parmap->getB("SAMPLER_BJORKEN_2D",false);
	Csampler::BJORKEN_YMAX=parmap->getD("SAMPLER_BJORKEN_YMAX",1.0);
	Csampler::USE_POLE_MASS=parmap->getB("SAMPLER_USE_POLE_MASS",false);
	Csampler::INCLUDE_BULK_VISCOSITY=parmap->getB("SAMPLER_INCLUDE_BULK_VISCOSITY",false);
	Csampler::INCLUDE_SHEAR_VISCOSITY=parmap->getB("SAMPLER_INCLUDE_SHEAR_VISCOSITY",false);
	int it,isigma;
	hyperlist.clear();
	sampler.resize(NTF+1);
	for(it=0;it<=NTF;it++){
		sampler[it].resize(NSIGMAF);
		for(isigma=0;isigma<NSIGMAF;isigma++){
			sampler[it][isigma]=nullptr;
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
	Csampler *samplerptr=nullptr,*sampler_findT=nullptr;
	//double Omega0Sum=0.0;
	partlist->nparts=0;
	list<Chyper *>::iterator it;

	for(it=hyperlist.begin();it!=hyperlist.end();it++){
		hyper=*it;
		if(hyper->firstcall){
			if(FINDT){
				if(sampler_findT==nullptr)
					sampler_findT=new Csampler(0.140,hyper->sigma);
				sampler_findT->GetTfMuNH(hyper);
				CALCMU=false;
			}
			samplerptr=ChooseSampler(hyper);
			hyper->sampler=samplerptr;
			if(samplerptr->FIRSTCALL){
				samplerptr->GetNHMu0();
				samplerptr->CalcDensitiesMu0();
				samplerptr->FIRSTCALL=false;
			}
			if(CALCMU)
				samplerptr->GetMuNH(hyper);
			hyper->firstcall=false;
			samplerptr->CalcNHadronsEpsilonP(hyper);
		}
		np=hyper->sampler->MakeParts(hyper);
		nparts+=np;
		//printf("np=%d,nparts=%d\n",np,nparts);
		//Omega0Sum+=hyper->dOmega[0];
	}
	if(sampler_findT!=nullptr)
		delete sampler_findT;
	NEVENTS+=1;
	return nparts;
}

Csampler* CmasterSampler::ChooseSampler(Chyper *hyper){
	double T,sigma,del;
	int it,isigma;
	T=hyper->T0;
	it=floorl(T/DELTF);
	del=T-it*DELTF;
	if(randy->ran()<del/DELTF)
		it+=1;
	if(it<0)
		it=0;
	if(it>=NTF){
		it=NTF-1;
		printf("WARNING in CmasterSampler::ChooseSampler\n");
		printf("your temperature, T=%g, is higher than TFmax=%g\n",T,TFmax);
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
	if(sampler[it][isigma]==nullptr){
		//printf("making Csampler object, DELTF=%g, it=%d, Tf-T0=%g\n",DELTF,it,(it+0.5)*DELTF-hyper->T0);
		sampler[it][isigma]=new Csampler((it+0.5)*DELTF,SIGMAFmin+(isigma+0.5)*DELSIGMAF);
	}
	return sampler[it][isigma];
}

void CmasterSampler::MakeDummyHyper(){
	Chyper *hyper=new Chyper();
	hyperlist.push_back(hyper);
}

void CmasterSampler::ReadHyper(){
	string filename;
	Chyper *elem;
	int ielement=0,alpha,beta;
	double u0,ux,uy,uz,utau,ueta;
	double t,x,y,z,tau,eta;
	double udotdOmega,udotdOmega_music;
	double dOmega0,dOmegaX,dOmegaY,dOmegaZ,dOmegaTau,dOmegaEta;
	double epsilonf;
	double Tdec,muB,muS,muC;
	double PIbulk __attribute__((unused)), Pdec __attribute__((unused));
	double qmu0,qmu1,qmu2,qmu3;
	double rhoB;
	double **pivisc=new double *[4];
	for(alpha=0;alpha<4;alpha++){
		pivisc[alpha]=new double[4];
		for(beta=0;beta<4;beta++)
			pivisc[alpha][beta]=0.0;
	}

	nelements=0;
	filename=parmap->getS("HYPER_INFO_FILE",string("../hydrodata/surface_2D.dat"));
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"rb");
	if (fptr==NULL) {
		fprintf(stderr,"Can't open hyper info file\n");
		printf("Error %d \n", errno);
	}

	double TotalVolume=0.0;

	while(!feof(fptr)){
		elem=new Chyper();
		// read from binary file
		float array[34];
		fread(array, sizeof(array),1,fptr);
		tau = array[0];
		x = array[1];
		y = array[2];
		eta = array[3];
		dOmegaTau = array[4];
		dOmegaX = array[5];
		dOmegaY = array[6];
		dOmegaEta = array[7];
		utau = array[8];
		ux = array[9];
		uy = array[10];
		ueta = array[11];
		utau=sqrt(1.0+ux*ux+uy*uy+ueta*ueta);
		
		const double u_milne_sqr = utau * utau - ux * ux - uy * uy - ueta * ueta;
		if (std::abs(u_milne_sqr - 1.0) > 1.e-6) {
			printf("Warning at reading from MUSIC output: "
				"u_Milne (u_eta multiplied by tau) = %9.6f %9.6f %9.6f %9.6f"
					", u^2 == 1 is not fulfilled with error %12.8f.\n",
			utau, ux, uy, ueta, std::abs(u_milne_sqr - 1.0));
		}
		udotdOmega_music = tau * (dOmegaTau * utau +
			dOmegaX * ux + dOmegaY * uy + dOmegaEta * ueta / tau);

		// Transforming from Milne to Cartesian
		const double ch_eta = std::cosh(eta), sh_eta = std::sinh(eta);
		t = tau * ch_eta;
		z = tau * sh_eta;
		u0 = utau * ch_eta + ueta * sh_eta;
		uz = utau * sh_eta + ueta * ch_eta;

		const double usqr = u0 * u0 - ux * ux - uy * uy - uz * uz;
		if (std::abs(usqr - 1.0) > 1.e-3) {
			printf("u*u should be 1, u*u = %12.8f\n", usqr);
		}

		dOmega0 = tau * ch_eta * dOmegaTau - sh_eta * dOmegaEta;
		dOmegaX = -tau * dOmegaX;
		dOmegaY = -tau * dOmegaY;
		dOmegaZ = tau * sh_eta * dOmegaTau - ch_eta * dOmegaEta;

		udotdOmega = dOmega0 * u0 - dOmegaX * ux - dOmegaY * uy - dOmegaZ * uz;
		if (std::abs(udotdOmega - udotdOmega_music) > 1.e-4) {
			printf("u^mu * dsigma_mu should be invariant: %12.9f == %12.9f\n",
			udotdOmega, udotdOmega_music);
		}
 
		epsilonf = array[12]*HBARC; //was labeled Edec--guessed this was epsilon
		Tdec = array[13]*HBARC;
		muB  = array[14]*HBARC/Tdec;
		muS  = array[15]*HBARC/Tdec;
		muC  = array[16]*HBARC/Tdec;
		Pdec = array[17]*Tdec - epsilonf;
		
		pivisc[0][0]=array[18]*HBARC;
		pivisc[0][1]=pivisc[1][0]=array[19]*HBARC;
		pivisc[0][2]=pivisc[2][0]=array[20]*HBARC;
		pivisc[0][3]=pivisc[3][0]=array[21]*HBARC; // /tau;
		pivisc[1][1]=array[22]*HBARC;
		pivisc[1][2]=pivisc[2][1]=array[23]*HBARC;
		pivisc[1][3]=pivisc[3][1]=array[24]*HBARC;  // /tau;
		pivisc[2][2]=array[25]*HBARC;
		pivisc[2][3]=pivisc[3][2]=array[26]*HBARC; // /tau;
		pivisc[3][3]=array[27]*HBARC; // /(tau*tau);
        TransformPiTotz(pivisc, ch_eta, sh_eta);
		
		//printf("reading, trace=%g =? 0\n",pivisc[0][0]-pivisc[1][1]-pivisc[2][2]-pivisc[3][3]);

		PIbulk = array[28]*HBARC;   // GeV/fm^3
		rhoB = array[29];  // 1/fm^3

		qmu0 = array[30];
		qmu1 = array[31];
		qmu2 = array[32];
		qmu3 = array[33];

		if(parmap->getB("SAMPLER_BJORKEN_2D",false)){
			double YMAX_ratio=parmap->getD("SAMPLER_BJORKEN_YMAX",1.0)/parmap->getD("HYDRO_BJORKEN_YMAX",1.0);
			dOmega0*=YMAX_ratio;
			dOmegaX*=YMAX_ratio;
			dOmegaY*=YMAX_ratio;
			dOmegaZ*=YMAX_ratio;
		}
		if(udotdOmega >= 0.0) {
			elem->tau=tau;
			elem->dOmega[0]=dOmega0; 
			elem->dOmega[1]=dOmegaX; 
			elem->dOmega[2]=dOmegaY; 
			elem->dOmega[3]=dOmegaZ;

			elem->udotdOmega=udotdOmega;

			elem->r[0]=t;
			elem->r[1]=x;
			elem->r[2]=y;
			elem->r[3]=z;

			elem->u[0]=u0;
			elem->u[1]=ux;
			elem->u[2]=uy;
			elem->u[3]=uz;
			
			GetPitilde(pivisc,elem->pitilde,elem->u);

			elem->muB=muB+0.5*muC;
			elem->muS=muS+0.5*muC;
			elem->muI=muC;
			
			elem->epsilon=epsilonf;
			elem->T0=Tdec;
			elem->rhoB=rhoB;
			elem->rhoS=0.0;
			elem->rhoI=0.0;

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
	//printf("Exiting ReadHyper() happily, TotalVolume=%lf, nelements=%d\n",TotalVolume,nelements);

    // avoid memory leak
    for (int i = 0; i < 4; i++)
        delete[] pivisc[i];
    delete[] pivisc;
}

void CmasterSampler::GetPitilde(double **pivisc,double **pitilde,FourVector &u){
	int alpha,beta;
	double picontract=0.0;
	double pivec[4]={0.0};
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pitilde[alpha][beta]=0.0;
			picontract+=pivisc[alpha][beta]*u[alpha]*u[beta];
			pivec[alpha]+=pivisc[alpha][beta]*u[beta];
		}
	}
	// Test Tr pivisc
#ifdef __TEST_PITILDE_TRACE__
	double trace=picontract/(u[0]*u[0]);
	for(alpha=1;alpha<4;alpha++){
		trace-=pivisc[alpha][alpha];
	}
	printf("-----------------------------\n");
	printf("before, trace=%g, pi_00=%g=?%g, u=(%g,%g,%g,%g)\n",trace,picontract/(u[0]*u[0]),pivisc[0][0],u[0],u[1],u[2],u[3]);
	printf("pivisc=\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%10.7f ",pivisc[alpha][beta]);
		}
		printf("\n");
	}
#endif
	//
	double x=u[0]*(1.0+u[0]);
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pitilde[alpha][beta]=pivisc[alpha][beta]-(u[alpha]*pivec[beta]/x)-(pivec[alpha]*u[beta]/x)+picontract*u[alpha]*u[beta]/(x*x);
		}
	}
#ifdef __TEST_PITILDE_TRACE__
	// Test Tr pitilde
	trace=0.0;
	for(alpha=1;alpha<4;alpha++){
		trace-=pitilde[alpha][alpha];
	}
	printf("after, trace=%g\n",trace);
	if(fabs(picontract)>0.1)
		Misc::Pause();
#endif
}

void CmasterSampler::TransformPiTotz(double **piMline, const double cosh_eta,
                                     const double sinh_eta) {
    double piCart[4][4];
    piCart[0][0] = (piMline[0][0]*cosh_eta*cosh_eta
                    + 2.*piMline[0][3]*cosh_eta*sinh_eta
                    + piMline[3][3]*sinh_eta*sinh_eta);
    piCart[0][1] = piMline[0][1]*cosh_eta + piMline[1][3]*sinh_eta;
    piCart[0][2] = piMline[0][2]*cosh_eta + piMline[2][3]*sinh_eta;
    piCart[0][3] = (piMline[0][0]*cosh_eta*sinh_eta
                    + piMline[0][3]*(cosh_eta*cosh_eta + sinh_eta*sinh_eta)
                    + piMline[3][3]*sinh_eta*cosh_eta);
    piCart[1][0] = piCart[0][1];
    piCart[1][1] = piMline[1][1];
    piCart[1][2] = piMline[1][2];
    piCart[1][3] = piMline[0][1]*sinh_eta + piMline[1][3]*cosh_eta;
    piCart[2][0] = piCart[0][2];
    piCart[2][1] = piCart[1][2];
    piCart[2][2] = piMline[2][2];
    piCart[2][3] = piMline[0][2]*sinh_eta + piMline[2][3]*cosh_eta;
    piCart[3][0] = piCart[0][3];
    piCart[3][1] = piCart[1][3];
    piCart[3][2] = piCart[2][3];
    piCart[3][3] = (piMline[0][0]*sinh_eta*sinh_eta
                    + 2.*piMline[0][3]*cosh_eta*sinh_eta
                    + piMline[3][3]*cosh_eta*cosh_eta);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            piMline[i][j] = piCart[i][j];
}
