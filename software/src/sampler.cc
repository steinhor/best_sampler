#include "pratt_sampler/sampler.h"
#include "pratt_sampler/sf.h"
#include "pratt_sampler/constants.h"


using namespace std;
using namespace pratt_sampler;

Crandy* Csampler::randy=NULL;
CresList *Csampler::reslist=NULL;
CparameterMap *Csampler::parmap=NULL;
CmasterSampler *Csampler::mastersampler=NULL;
bool Csampler::CALCMU=false;
int Csampler::n_bose_corr=1;
bool Csampler::bose_corr=false;
bool Csampler::BJORKEN_2D=false;
double Csampler::BJORKEN_YMAX=1.0;
bool Csampler::USE_POLE_MASS=false;

// Constructor
Csampler::Csampler(double Tfset,double sigmafset){
	FIRSTCALL=true;
	Tf=Tfset;
	sigmaf=sigmafset;
	if(!bose_corr)
		n_bose_corr=1;
	int nres=reslist->massmap.size();
	if(reslist->GetResInfoPtr(22)->pid==22)
		nres-=1;
	density0i.resize(nres);
	epsilon0i.resize(nres);
	P0i.resize(nres);
	lambda0i.resize(nres);
	sfdens0imap.resize(nres);
	if(bose_corr){
		pibose_P0.resize(n_bose_corr+1);
		pibose_epsilon0.resize(n_bose_corr+1);
		pibose_dedt0.resize(n_bose_corr+1);
		pibose_dens0.resize(n_bose_corr+1);
		pibose_lambda0.resize(n_bose_corr+1);
	}
	bose_test_off=false;
	bose_test=false;
	viscous_test=false;
}

// Destructor
Csampler::~Csampler(){
}

// Calculates array of densities for each resonance with mu=0, to be used in sampling
void Csampler::CalcDensitiesMu0(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,dummy;
	int ires;
	nhadrons0=P0=epsilon0=0.0;

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			GetDensPMax(resinfo,densi,epsiloni,Pi,dedti);
			ires=resinfo->ires;
			density0i[ires]=densi;
			epsilon0i[ires]=epsiloni;
			P0i[ires]=Pi;
			nhadrons0+=density0i[ires];
			epsilon0+=epsiloni;
			P0+=Pi;
			if(resinfo->decay)
				CalcSFDensMap(resinfo,Tf,sfdens0imap[ires]);
		}
	}
	if(bose_corr){
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			EOS::freegascalc_onespecies(Tf/nbose,0.138,pibose_epsilon0[nbose],pibose_P0[nbose],pibose_dens0[nbose],
			dummy,pibose_dedt0[nbose]);
			nhadrons0+=3.0*pibose_dens0[nbose];
			epsilon0+=3.0*pibose_epsilon0[nbose];
			P0+=3.0*pibose_P0[nbose];
		}
	}
}

// Calculates factors (depend only on T) used for Newton's method to get muB, muI, muS from rhoB, rhoI, rhoS
void Csampler::GetNHMu0(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double Pi,epsiloni,densi,dedti,dummy;
	int B,S,II,Q;

	nh0_b0i0s0=nh0_b0i2s0=nh0_b0i1s1=0.0;
	nh0_b1i0s1=nh0_b1i0s3=nh0_b1i1s0=nh0_b1i1s2=nh0_b1i2s1=nh0_b1i3s0=0.0;
	nh0_b2i0s0=0.0;

	eh0_b0i0s0=eh0_b0i2s0=eh0_b0i1s1=0.0;
	eh0_b1i0s1=eh0_b1i0s3=eh0_b1i1s0=eh0_b1i1s2=eh0_b1i2s1=eh0_b1i3s0=0.0;
	eh0_b2i0s0=0.0;

	dedth0_b0i0s0=dedth0_b0i2s0=dedth0_b0i1s1=0.0;
	dedth0_b1i0s1=dedth0_b1i0s3=dedth0_b1i1s0=dedth0_b1i1s2=dedth0_b1i2s1=dedth0_b1i3s0=0.0;
	dedth0_b2i0s0=0.0;

	if(bose_corr){
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			EOS::freegascalc_onespecies(Tf/nbose,0.138,pibose_epsilon0[nbose],pibose_P0[nbose],pibose_dens0[nbose],
			dummy,pibose_dedt0[nbose]);
		}
	}

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			GetDensPMax(resinfo,densi,epsiloni,Pi,dedti);

			B=resinfo->baryon;
			S=resinfo->strange;
			Q=resinfo->charge;
			II=2*Q-B-S; // II is 2*I3
			B=abs(B);
			S=abs(S);
			II=abs(II);

			if(B==0 && II==0 && S==0){
				nh0_b0i0s0+=densi;
				eh0_b0i0s0+=epsiloni;
				dedth0_b0i0s0+=dedti;
			}
			else if(B==0 && II==1 && S==1){
				nh0_b0i1s1+=densi;
				eh0_b0i1s1+=epsiloni;
				dedth0_b0i1s1+=dedti;
			}
			else if(B==0 && II==2 && S==0){
				nh0_b0i2s0+=densi;
				eh0_b0i2s0+=epsiloni;
				dedth0_b0i2s0+=dedti;
			}
			else if(B==1 && II==0 && S==1){
				nh0_b1i0s1+=densi;
				eh0_b1i0s1+=epsiloni;
				dedth0_b1i0s1+=dedti;
			}
			else if(B==1 && II==0 && S==3){
				nh0_b1i0s3+=densi;
				eh0_b1i0s3+=epsiloni;
				dedth0_b1i0s3+=dedti;
			}
			else if(B==1 && II==1 && S==0){
				nh0_b1i1s0+=densi;
				eh0_b1i1s0+=epsiloni;
				dedth0_b1i1s0+=dedti;
			}
			else if(B==1 && II==1 && S==2){
				nh0_b1i1s2+=densi;
				eh0_b1i1s2+=epsiloni;
				dedth0_b1i1s2+=dedti;
			}
			else if(B==1 && II==2 && S==1){
				nh0_b1i2s1+=densi;
				eh0_b1i2s1+=epsiloni;
				dedth0_b1i2s1+=dedti;
			}
			else if(B==1 && II==3 && S==0){
				nh0_b1i3s0+=densi;
				eh0_b1i3s0+=epsiloni;
				dedth0_b1i3s0+=dedti;
			}
			else if(B==2 && II==0 && S==0){ //deuteron, not currently used in calculations
				nh0_b2i0s0+=densi;
				eh0_b2i0s0+=epsiloni;
				dedth0_b2i0s0+=dedti;
			}
			else{
				printf("NO BIS Match!!!! B=%d, II=%d, S=%d\n",B,II,S);
				exit(1);
			}
		}
	}
}

// These routines are used for calculating the factor used for viscous corrections, lambda
// Factor used in CalcLambda
double Csampler::GetDIppForLambda(CresInfo *resinfo,double TT){
	double m,degen,z,J,nfact,sign,temp;
	double	dIpp=0.0;
	const int nmax=70;
	double G[nmax+5];
	int n;
	m=resinfo->mass;
	degen=resinfo->degen;
	z=m/TT;
	G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
	for(int j=1;j<nmax+5;j++){
		n=5-2*j;
		if(n!=-1)	G[j]=(-exp(-z)/n)+(G[j-1]*z*z-z*exp(-z))/((n+1.0)*n);
		else G[j]=gsl_sf_gamma_inc(-1,z)*z;
	}
	J=0.0;
	nfact=1.0;
	sign=1.0;
	for(n=0;n<nmax;n+=1){
		if(n>0) sign=-1.0;
		J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
		nfact=nfact*0.5/(n+1.0);
		if(n>0) nfact*=(2.0*n-1.0);
	}
	temp=degen*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
	dIpp+=temp;
	dIpp=dIpp/(60.0*M_PI*M_PI*HBARC*HBARC*HBARC);
	return dIpp;
}
// for Mu=0
void Csampler::CalcLambdaMu0(){
	double Ipp=0.0,dIpp;
	int nbose;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->pid!=22){
			dIpp=GetDIppForLambda(resinfo,Tf);
			lambda0i[resinfo->ires]=dIpp;
			Ipp+=dIpp;
		}
	}
	if(bose_corr){
		for(nbose=2;nbose<=n_bose_corr;nbose++){
			dIpp=GetDIppForLambda(resinfo,Tf/double(nbose));
			pibose_lambda0[nbose]=dIpp;
			Ipp+=dIpp;
		}
	}
	lambda0=Ipp;
}
// for Mu!=0
double Csampler::CalcLambdaF(Chyper *hyper){
	return CalcLambdaF(hyper->muB,hyper->muI,hyper->muS,hyper->P);
}
double Csampler::CalcLambdaF(double muB,double muI,double muS,double Pf){
	int ires=0,nbose;
	double Ipp=0.0,dIpp,xx,lambdaf;
	double lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	if(mastersampler->SETMU0){
		return lambda0;
	}
	else{

		for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
			resinfo=rpos->second;
			if(resinfo->pid!=22){
				I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
				mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
				dIpp=0.0;
				Ipp+=lambda0i[ires]*exp(mutot);
				ires+=1;
			}
		}
		if(bose_corr){
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				xx=exp(muI);
				Ipp+=pibose_lambda0[nbose]*(pow(xx,nbose)+pow(xx,-nbose)+1.0);
			}
		}
		if(mastersampler->SETMU0)
			lambdafact=2.0*P0-4.0*Ipp;
		else
			lambdafact=2.0*Pf-4.0*Ipp;
		lambdaf=lambdafact;
		return lambdaf;
	}
}

// Calculates muB, muI, muS from rhoB, rhoI, rhoS, also calculates nhadronsf (factors above must be calculated first)
void Csampler::GetMuNH(Chyper *hyper){
	GetMuNH(hyper->rhoB,hyper->rhoI,hyper->rhoS,hyper->muB,hyper->muI,hyper->muS);
}
void Csampler::GetMuNH(double rhoBtarget,double rhoItarget,double rhoStarget,double &muB,double &muI,double &muS){
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd mu(3),dmu(3),drho(3);

	//3D Newton's Method
	// Here rhoI refers to rho_u-rho_d = 2*I3 and mu[1]=muI/2
	double rhoB,rhoS,rhoI;

	double drhoB_dmuB,drhoB_dmuS,drhoB_dmuI;
	double drhoS_dmuB,drhoS_dmuS,drhoS_dmuI;
	double drhoI_dmuB,drhoI_dmuS,drhoI_dmuI;

	double xB,xI,xS,xxB,xxI,xxS,dmumag;
	int ntries=0;

	mu[0]=muB;
	mu[1]=0.5*muI;
	mu[2]=muS;
	do{
		ntries+=1;
		xB=exp(mu[0]);
		xI=exp(mu[1]);
		xS=exp(mu[2]);
		xxB=1.0/xB;
		xxI=1.0/xI;
		xxS=1.0/xS;

		rhoB=0.5*nh0_b1i0s1*(xB*xxS-xxB*xS)
			+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
				+0.25*nh0_b1i1s0*(xB-xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
							+0.25*nh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI);
		drhoB_dmuB=0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
		drhoB_dmuI=0.25*nh0_b1i1s0*(xB-xxB)*(xI-xxI)
			+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI-xxI)
				+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(2*xI*xI-2*xxI*xxI)
					+0.25*nh0_b1i3s0*(xB-xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
		drhoB_dmuS=0.5*nh0_b1i0s1*(-xB*xxS-xxB*xS)
			+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS-3*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS-2*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(-xB*xxS-xxB*xS)*(xI*xI+xxI*xxI);

		rhoI=0.5*nh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
			+0.25*nh0_b0i1s1*(xI-xxI)*(xS+xxS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI-xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
		drhoI_dmuB=drhoB_dmuI;
		drhoI_dmuI=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
					+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
							+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);
		if (bose_corr) {
			for(int nbose=2;nbose<=n_bose_corr;nbose++){
				rhoI+=pibose_dens0[nbose]*2.0*sinh(2.0*nbose*mu[1]);
				drhoI_dmuI+=pibose_dens0[nbose]*nbose*4.0*cosh(2.0*nbose*mu[1]);
			}
		}
		drhoI_dmuS=0.25*nh0_b0i1s1*(xI-xxI)*(xS-xxS)
			+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI-xxI)
				+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI);

		rhoS=0.25*nh0_b0i1s1*(xI+xxI)*(xS-xxS)
			+0.5*nh0_b1i0s1*(-xB*xxS+xxB*xS)
				+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
					+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);
		drhoS_dmuB=drhoB_dmuS;
		drhoS_dmuI=drhoI_dmuS;
		drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
					+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
						+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

		drho[0]=rhoBtarget-rhoB;
		drho[1]=rhoItarget-rhoI;
		drho[2]=rhoStarget-rhoS;

		A(0,0)=drhoB_dmuB;
		A(0,1)=A(1,0)=drhoB_dmuI;
		A(0,2)=A(2,0)=drhoB_dmuS;
		A(1,1)=drhoI_dmuI;
		A(1,2)=A(2,1)=drhoI_dmuS;
		A(2,2)=drhoS_dmuS;
		//cout << A << endl;
		dmu=A.colPivHouseholderQr().solve(drho);
		dmumag=sqrt(dmu[0]*dmu[0]+dmu[1]*dmu[1]+dmu[2]*dmu[2]);
		if(dmumag>2.0){
			dmu[0]=2.0*dmu[0]/dmumag;
			dmu[1]=2.0*dmu[1]/dmumag;
			dmu[2]=2.0*dmu[2]/dmumag;
		}
		mu[0]+=dmu[0]; mu[1]+=dmu[1]; mu[2]+=dmu[2];

	}while(fabs(drho[0])>1.0E-8 || fabs(drho[1])>1.0E-8 || fabs(drho[2])>1.0E-8);
	muB=mu[0];
	muI=2.0*mu[1];
	muS=mu[2];
	//printf("rho=(%g,%g,%g) =? (%g,%g,%g)\n",rhoB,rhoI,rhoS,rhoBtarget,rhoItarget,rhoStarget);
}

// Same as above, but also calculates T, also uses epsilon (uses GetEpsilonRhoDerivatives to get factors )
void Csampler::GetTfMuNH(Chyper *hyper){
	GetTfMuNH(hyper->epsilon,hyper->rhoB,hyper->rhoI,hyper->rhoS,hyper->muB,hyper->muI,hyper->muS);
	hyper->T=Tf;
}
void Csampler::GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget,double &muB,double &muI,double &muS){
	//4D Newton's Method
	// Here rhoI refers to rho_u-rho_d = 2*I3 and mu[1]=muI/2
	double epsilon=0,rhoB=0,rhoS=0,rhoI=0;
	Eigen::MatrixXd A(4,4);
	Eigen::VectorXd dmu(4),drho(4);
	double cmb,smb;
	int ntries=0;
	do{
		ntries+=1;
		if(ntries>30000){
			printf("FAILURE, ntries=%d\n",ntries);
			exit(1);
		}
		smb=sinh(muB);
		cmb=cosh(muB);

		GetNHMu0();
		GetEpsilonRhoDerivatives(muB,muI,muS,epsilon,rhoB,rhoI,rhoS,A);
		for(int i=0;i<4;i++){
			A(i,1)=A(i,1)/cmb;
		}
		drho[0]=epsilontarget-epsilon;
		drho[1]=rhoBtarget-rhoB;
		drho[2]=rhoItarget-rhoI;
		drho[3]=rhoStarget-rhoS;
		dmu=A.colPivHouseholderQr().solve(drho);
		Tf+=dmu[0];
		smb+=dmu[1];
		muB=asinh(smb);
		muI+=dmu[2];
		muS+=dmu[3];
	}while(fabs(drho[0])>1.0E-4 || fabs(drho[1])>1.0E-6 || fabs(drho[2])>1.0E-6 || fabs(drho[3])>1.0E-6);
	//printf("Tf=%g, Mu=(%g,%g,%g), epsilon=%g=?%g, rho=(%g,%g,%g) =? (%g,%g,%g)\n",
	//Tf,muB,muI,muS,epsilon,epsilontarget,rhoB,rhoI,rhoS,rhoBtarget,rhoItarget,rhoStarget);
}
void Csampler::GetEpsilonRhoDerivatives(double muB,double muI,double muS,double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A){
	double xB,xI,xS,xxB,xxI,xxS;

	double drhoB_dT,drhoB_dmuB,drhoB_dmuS,drhoB_dmuI;
	double drhoS_dT,drhoS_dmuB,drhoS_dmuS,drhoS_dmuI;
	double drhoI_dT,drhoI_dmuB,drhoI_dmuS,drhoI_dmuI;
	double de_dT,de_dmuB,de_dmuI,de_dmuS;

	GetNHMu0();

	xB=exp(muB);
	xI=exp(0.5*muI);
	xS=exp(muS);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;

	epsilon=eh0_b0i0s0+0.5*eh0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*eh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*eh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*eh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*eh0_b2i0s0*(xB*xB+xxB*xxB);

	de_dT=dedth0_b0i0s0+0.5*dedth0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*dedth0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*dedth0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*dedth0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*dedth0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*dedth0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*dedth0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*dedth0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
									+0.5*dedth0_b2i0s0*(xB*xB+xxB*xxB);


	de_dmuB=0.5*eh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*eh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*eh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+eh0_b2i0s0*(xB*xB-xxB*xxB);

	de_dmuI=0.5*eh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
		+0.25*eh0_b0i1s1*(xI-xxI)*(xS+xxS)
			+0.25*eh0_b1i1s0*(xB+xxB)*(xI-xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
					+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
						+0.25*eh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);

	de_dmuS=0.25*eh0_b0i1s1*(xI+xxI)*(xS-xxS)
		+0.5*eh0_b1i0s1*(-xB*xxS+xxB*xS)
			+0.5*eh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
				+0.25*eh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	rhoB=0.5*nh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*nh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*nh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+nh0_b2i0s0*(xB*xB-xxB*xxB);

	drhoB_dT=de_dmuB/(Tf*Tf);

	drhoB_dmuB=0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
		+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
							+2.0*nh0_b2i0s0*(xB*xB-xxB*xxB);

	drhoB_dmuI=0.25*nh0_b1i1s0*(xB-xxB)*(xI-xxI)
		+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI-xxI)
			+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(2*xI*xI-2*xxI*xxI)
				+0.25*nh0_b1i3s0*(xB-xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);

	drhoB_dmuS=0.5*nh0_b1i0s1*(-xB*xxS-xxB*xS)
		+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS-3*xxB*xS*xS*xS)
			+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS-2*xxB*xS*xS)*(xI+xxI)
				+0.25*nh0_b1i2s1*(-xB*xxS-xxB*xS)*(xI*xI+xxI*xxI);

	rhoI=0.5*nh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
		+0.25*nh0_b0i1s1*(xI-xxI)*(xS+xxS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI-xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);

	drhoI_dT=de_dmuI/(Tf*Tf);

	drhoI_dmuB=drhoB_dmuI;

	drhoI_dmuI=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
						+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);

	drhoI_dmuS=0.25*nh0_b0i1s1*(xI-xxI)*(xS-xxS)
		+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI-xxI)
			+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI);

	rhoS=0.25*nh0_b0i1s1*(xI+xxI)*(xS-xxS)
		+0.5*nh0_b1i0s1*(-xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	drhoS_dT=de_dmuS/(Tf*Tf);
	drhoS_dmuB=drhoB_dmuS;
	drhoS_dmuI=drhoI_dmuS;
	drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
		+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	if (bose_corr) {
		double ch,sh;
		for(int nbose=2;nbose<=n_bose_corr;nbose++){
			ch=cosh(nbose*muI);
			sh=sinh(nbose*muI);
			epsilon+=pibose_epsilon0[nbose]*(2.0*ch+1.0);
			de_dT+=pibose_dedt0[nbose]*(2.0*ch+1.0);
			de_dmuI+=pibose_epsilon0[nbose]*4.0*nbose*sh;
			rhoI+=pibose_dens0[nbose]*2.0*sh;
			drhoI_dT+=pibose_epsilon0[nbose]*4.0*nbose*sh/(Tf*Tf);
			drhoI_dmuI+=pibose_dens0[nbose]*nbose*4.0*ch;
		}
	}

	A(0,0)=de_dT;
	A(0,1)=de_dmuB;
	A(0,2)=0.5*de_dmuI;
	A(0,3)=de_dmuS;

	A(1,0)=drhoB_dT;
	A(1,1)=drhoB_dmuB;
	A(1,2)=0.5*drhoB_dmuI;
	A(1,3)=drhoB_dmuS;

	A(2,0)=drhoI_dT;
	A(2,1)=drhoI_dmuB;
	A(2,2)=0.5*drhoI_dmuI;
	A(2,3)=drhoI_dmuS;

	A(3,0)=drhoS_dT;
	A(3,1)=drhoS_dmuB;
	A(3,2)=0.5*drhoS_dmuI;
	A(3,3)=drhoS_dmuS;

	//printf("Calculating rho,epsilon,A: T=%g, epsilon=%g, de_dT=%g\n",Tf,epsilon,de_dT);
	//printf("mu=(%g,%g,%g), rho=(%g,%g,%g)\n",muB,muI,muS,rhoB,rhoI,rhoS);
	//cout << A << endl;
	//exit(1);
}

// For sampling, generates mass of resonances (uses GetDensPMax below)
double Csampler::GenerateThermalMass(CresInfo *resinfo){
	map<double,double>::iterator it1,it2;
	double E,E1,E2;
	int ires=resinfo->ires;
	double p1,p2,netprob=randy->ran();
	if(!resinfo->decay)
		E=resinfo->mass;
	else{
		it1=sfdens0imap[ires].lower_bound(netprob);
		if(it1==sfdens0imap[ires].end()){
			printf("it1 already at end of map\n");
			exit(1);
		}
		--it1;
		it2=it1; it2++;
		p1=it1->first;
		p2=it2->first;
		E1=it1->second;
		E2=it2->second;
		E=((netprob-p1)*E2+(p2-netprob)*E1)/(p2-p1);
	}
	if(E<resinfo->minmass || E1>E || E2<E){
		printf("something odd, E=%g, (E1,E2)=(%g,%g), (p1,p2)=(%g,%g), minmass=%g, netprob=%g\n",E,E1, E2,p1,p2,resinfo->minmass,netprob);
		resinfo->PrintBranchInfo();
		resinfo->PrintSpectralFunction();
		map<double,double>::iterator it;
		for(it=sfdens0imap[ires].begin();it!=sfdens0imap[ires].end();++it){
			printf("Emap=%g, netdens=%g\n",it->second,it->first);
		}
		exit(1);
	}
	return E; //returns a random mass proportional to dens*SF'
}

// Calculates density, pressure, energy density, de/dt (for MC mass choice) for individual resonances
void Csampler::GetDensPMax(CresInfo *resinfo,double &densi,double &epsiloni,double &Pi,double &dedti){
	double degen,width,m,minmass,sigma2;
	CmeanField *mf=mastersampler->meanfield;
	bool decay=resinfo->decay;
	//decay=false;
	degen=resinfo->degen;
	m=mf->GetMass(resinfo,sigmaf);
	width=resinfo->width;
	minmass=resinfo->minmass;

	if((minmass>-0.001) && (width>0.000000001) && decay && !USE_POLE_MASS){
		EOS::freegascalc_onespecies_finitewidth(resinfo,Tf,epsiloni,Pi,densi,sigma2,dedti);
	}
	else{
		EOS::freegascalc_onespecies(Tf,m,epsiloni,Pi,densi,sigma2,dedti);
	}
	epsiloni*=degen;
	Pi*=degen;
	dedti*=degen;
	densi*=degen;
}

// gets nhadrons, epsilon and P
void Csampler::GetNHadronsEpsilonPF(Chyper *hyper){
	GetNHadronsEpsilonPF(hyper->muB,hyper->muI,hyper->muS,hyper->nhadrons,hyper->epsilon,hyper->P);
}

void Csampler::GetNHadronsEpsilonPF(double muB,double muI,double muS,double &nhadronsf,double &epsilonf,double &Pf){
	double xB,xI,xS,xxB,xxI,xxS,xbose;
	int nbose;
	if(mastersampler->SETMU0){
		nhadronsf=nhadrons0;
		epsilonf=epsilon0;
		Pf=P0;
	}
	else{
		xB=exp(muB);
		xI=exp(muI);
		xS=exp(muS);
		xxB=1.0/xB;
		xxI=1.0/xI;
		xxS=1.0/xS;
		xbose=exp(muI);
		nhadronsf=nh0_b0i0s0+0.5*nh0_b0i2s0*(xI*xI+xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
										+0.5*nh0_b2i0s0*(xB*xB+xxB*xxB);
		Pf=nhadronsf*Tf;
		epsilonf=eh0_b0i0s0+0.5*eh0_b0i2s0*(xI*xI+xxI*xxI)
			+0.25*eh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*eh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*eh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI)
										+0.5*eh0_b2i0s0*(xB*xB+xxB*xxB);

		if (bose_corr){
			for(nbose=2;nbose<=n_bose_corr;nbose++){
				xbose=(xbose+1.0+1.0/xbose);
				nhadronsf+=pibose_dens0[nbose]*xbose;
				epsilonf+=pibose_epsilon0[nbose]*xbose;
				Pf+=pibose_P0[nbose]*(pow(xbose,nbose)+pow(xbose,-nbose)+1.0);
			}
		}
	}
}

void Csampler::CalcSFDensMap(CresInfo *resinfo,double T,map<double,double> &sfdensmap){
	int iE,nE;
	double z,E,E2,sfnorm,netprob,k1,k0;
	vector<double> dens;
	nE=resinfo->SpectVec.size();
	dens.resize(nE);
	sfnorm=0.0;
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		z=E/T;
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		dens[iE]=E*E*k0+2.0*E*T*k1;
		sfnorm+=resinfo->SpectVec[iE]*dens[iE];
	}
	netprob=0.0;
	sfdensmap.insert(pair<double,double>(netprob,resinfo->minmass));
	for(iE=0;iE<nE;iE++){
		E=resinfo->SpectEVec[iE];
		netprob+=resinfo->SpectVec[iE]*dens[iE]/sfnorm;
		if(iE!=nE-1)
			E2=0.5*(resinfo->SpectEVec[iE]+resinfo->SpectEVec[iE+1]);
		else
			E2=2.0*resinfo->SpectEVec[iE]-resinfo->SpectEVec[iE-1];
		sfdensmap.insert(pair<double,double>(netprob,E2));
	}
	dens.clear();
}
