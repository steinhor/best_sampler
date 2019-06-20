#include "sampler.h"
using namespace std;
Crandy* Csampler::randy=NULL;
CresList *Csampler::reslist=NULL;
CparameterMap *Csampler::parmap=NULL;
CmasterSampler *Csampler::mastersampler=NULL;

Csampler::Csampler(double Tfset,double sigmafset){
    Tf=Tfset;
	sigmaf=sigmafset;
	muB=muS=muI=0.0;
	int nres=reslist->massmap.size();
	if(reslist->GetResInfoPtr(22)->code==22)
		nres-=1;
	densityf0.resize(nres);
	densityf.resize(nres);
	epsilonf0i.resize(nres);
	Pf0i.resize(nres);
	lambdaf0i.resize(nres);
	maxweight.resize(nres);
	CalcDensitiesF0();
	CalcLambdaF0();
	GetNH0();
}

Csampler::~Csampler(){
}

void Csampler::CalcLambda(){
	int iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign;
	double dIpptest=0.0,dp=4.0,p,e,lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=resinfo->spin;
			z=m/Tf;
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-muS*resinfo->strange);
			mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(int i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
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
			dIpp=degen*exp(mutot)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			Ipp+=dIpp;

		}
	}
	if(mastersampler->SETMU0)
		lambdafact=2.0*Pf0-4.0*Ipp;
	else
		lambdafact=2.0*Pf-4.0*Ipp;
	lambdaf=lambdafact;
}

void Csampler::CalcLambdaF0(){
	int iQ,n,i,ires;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign;
	double dIpptest=0.0,dp=4.0,p,e,lambdafact,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	ires=0;
	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=resinfo->spin;
			z=m/Tf;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(int i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
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
			dIpp=degen*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			lambdaf0i[resinfo->ires]=dIpp;
			Ipp+=dIpp;
			ires+=1;
		}
	}
	if(mastersampler->SETMU0)
		lambdafact=2.0*Pf0-4.0*Ipp;
	else
		lambdafact=2.0*Pf-4.0*Ipp;
	lambdaf0=lambdafact;
	lambdaf=lambdaf0;
}

void Csampler::CalcLambdaF(){
	int iQ,n,ires;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign;
	double dIpptest=0.0,dp=4.0,p,e,lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	ires=0;
	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-muS*resinfo->strange);
			mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
			Ipp+=lambdaf0i[ires]*exp(mutot);
			ires+=1;
		}
	}
	if(mastersampler->SETMU0)
		lambdafact=2.0*Pf0-4.0*Ipp;
	else
		lambdafact=2.0*Pf-4.0*Ipp;
	lambdaf=lambdafact;
}

void Csampler::GetNH0(){
    CresInfo *resinfo;
    CresMassMap::iterator rpos;
    double m,degen;
    double Pi,epsiloni,densi,sigma2i,dedti,si,maxweighti;
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
    nhadronsf0=0.0;

    for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
        resinfo=rpos->second;
        if(resinfo->code!=22){
            GetDensPMaxWeight(resinfo,0.0,densi,epsiloni,Pi,dedti,maxweighti);
            nhadronsf0+=densi;

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
            else if(B==2 && II==0 && S==0){ //deuteron
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

void Csampler::GetMuNH(Chyper *hyper){
	// Here rhoI refers to rho_u-rho_d = 2*I3 and mu[1]=muI/2
	double rhoBtarget=hyper->rhoB;
	double rhoItarget=2.0*hyper->rhoI;
	double rhoStarget=0.0;
	double rhoB,rhoS,rhoI;
	double drhoB_dmuB,drhoB_dmuS,drhoB_dmuI;
	double drhoS_dmuB,drhoS_dmuS,drhoS_dmuI;
	double drhoI_dmuB,drhoI_dmuS,drhoI_dmuI;
	double xB,xI,xS,xxB,xxI,xxS;
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd mu(3),dmu(3),drho(3);
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
		//printf("drho=(%g,%g,%g)\n",drho[0],drho[1],drho[2]);
		A(0,0)=drhoB_dmuB;
		A(0,1)=A(1,0)=drhoB_dmuI;
		A(0,2)=A(2,0)=drhoB_dmuS;
		A(1,1)=drhoI_dmuI;
		A(1,2)=A(2,1)=drhoI_dmuS;
		A(2,2)=drhoS_dmuS;
		dmu=A.colPivHouseholderQr().solve(drho);
		//printf("ntries=%d: rhoB=%g, rhoBtarget=%g, rhoI=%g, rhoItarget=%g, rhoS=%g\n",ntries,rhoB, rhoBtarget,rhoI,rhoItarget,rhoS);
		mu[0]+=dmu[0]; mu[1]+=dmu[1]; mu[2]+=dmu[2];

	}while(fabs(drho[0])>1.0E-8 || fabs(drho[1])>1.0E-8 || fabs(drho[2])>1.0E-8);

	xB=exp(mu[0]);
	xI=exp(mu[1]);
	xS=exp(mu[2]);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;
	nhadronsf=nh0_b0i0s0+0.5*nh0_b0i2s0*(xI*xI*+xxI*xxI)
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
	epsilonf=eh0_b0i0s0+0.5*eh0_b0i2s0*(xI*xI*+xxI*xxI)
		+0.25*eh0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*eh0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*eh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
	hyper->muB=mu[0];
	hyper->muI=2.0*mu[1];
	hyper->muS=mu[2];
	//printf("mu=(%g,%g,%g), rho=(%g,%g,%g) =? (%g,%g,0)\n",muB,muI,muS,rhoB,rhoI,rhoS,rhoBtarget,rhoItarget);
	CalcDensitiesF();
}

void Csampler::GetTfMuNH(Chyper *hyper){
	Tf=hyper->T;
	muB=hyper->muB;
	muI=hyper->muI;
	muS=hyper->muS;
	GetTfMuNH(hyper->epsilon,hyper->rhoB,hyper->rhoI,hyper->rhoS);
	hyper->T=Tf;
	hyper->muB=muB;
	hyper->muI=muI;
	hyper->muS=muS;
	printf("TEST: nhadronsf=%g\n",nhadronsf);
	hyper->nhadrons=nhadronsf;
}

void Csampler::GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget){
	// Here rhoI refers to rho_u-rho_d = 2*I3 and mu[1]=muI/2
	double epsilon,rhoB,rhoS,rhoI;
	double xB,xI,xS,xxB,xxI,xxS;
	double alpha=1.0;
	int i;
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
		GetEpsilonRhoDerivatives(epsilon,rhoB,rhoI,rhoS,A);
		for(i=0;i<4;i++){
			A(i,1)=A(i,1)/cmb;
			//A(i,0)=A(i,0)/(alpha*pow(Tf,alpha-1.0));
		}
		drho[0]=epsilontarget-epsilon;
		drho[1]=rhoBtarget-rhoB;
		drho[2]=rhoItarget-rhoI;
		drho[3]=rhoStarget-rhoS;
		dmu=A.colPivHouseholderQr().solve(drho);
		//Tf=pow(pow(Tf,alpha)+dmu[0],1.0/alpha);
		Tf+=dmu[0];
		smb+=dmu[1];
		muB=asinh(smb);
		muI+=dmu[2];
		muS+=dmu[3];

	}while(fabs(drho[0])>1.0E-3 || fabs(drho[1])>1.0E-5 || fabs(drho[2])>1.0E-5 || fabs(drho[3])>1.0E-5);

	xB=exp(muB);
	xI=exp(0.5*muI);
	xS=exp(muS);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;
	//xB=xI=xS=xxB=xxI=xxS=1.0;
	nhadronsf=nh0_b0i0s0
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*nh0_b0i2s0*(xI*xI+xxI*xxI)
					+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);

	CalcDensitiesF();
}

void Csampler::GetEpsilonRhoDerivatives(double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A){
	double xB,xI,xS,xxB,xxI,xxS;
	double drhoB_dt,drhoB_dmuB,drhoB_dmuS,drhoB_dmuI;
	double drhoS_dt,drhoS_dmuB,drhoS_dmuS,drhoS_dmuI;
	double drhoI_dt,drhoI_dmuB,drhoI_dmuS,drhoI_dmuI;
	double de_dT,de_dmuB,de_dmuI,de_dmuS;
	GetNH0();
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
								+0.25*eh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);

	de_dT=dedth0_b0i0s0+0.5*dedth0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*dedth0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*dedth0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*dedth0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*dedth0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*dedth0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*dedth0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*dedth0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);

	de_dmuB=0.5*eh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*eh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*eh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*eh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI);
	de_dmuI=0.5*eh0_b0i2s0*(2*xI*xI-2*xxI*xxI)
		+0.25*eh0_b0i1s1*(xI-xxI)*(xS+xxS)
			+0.25*eh0_b1i1s0*(xB+xxB)*(xI-xxI)
				+0.25*eh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI-xxI)
					+0.25*eh0_b1i2s1*(xB*xxS+xxB*xS)*(2*xI*xI-2*xxI*xxI)
						+0.25*eh0_b1i3s0*(xB+xxB)*(3*xI*xI*xI-3*xxI*xxI*xxI);
	de_dmuS=0.25*eh0_b0i1s1*(xI+xxI)*(xS-xxS)
		+0.5*eh0_b1i0s1*(-xB*xxS-xxB*xS)
			+0.5*eh0_b1i0s3*(-3*xB*xxS*xxS*xxS+3*xxB*xS*xS*xS)
				+0.25*eh0_b1i1s2*(-2*xB*xxS*xxS+2*xxB*xS*xS)*(xI+xxI)
					+0.25*eh0_b1i2s1*(-xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	rhoB=0.5*nh0_b1i0s1*(xB*xxS-xxB*xS)
		+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS-xxB*xS*xS*xS)
			+0.25*nh0_b1i1s0*(xB-xxB)*(xI+xxI)
				+0.25*nh0_b1i1s2*(xB*xxS*xxS-xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS-xxB*xS)*(xI*xI+xxI*xxI)
						+0.25*nh0_b1i3s0*(xB-xxB)*(xI*xI*xI+xxI*xxI*xxI);
	drhoB_dt=de_dmuB/(Tf*Tf);
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
	drhoI_dt=de_dmuI/(Tf*Tf);
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

	drhoS_dt=de_dmuS/(Tf*Tf);
	drhoS_dmuB=drhoB_dmuS;
	drhoS_dmuI=drhoI_dmuS;
	drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
		+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

	//printf("drho=(%g,%g,%g)\n",drho[0],drho[1],drho[2]);
	A(0,0)=de_dT;
	A(0,1)=de_dmuB;
	A(0,2)=0.5*de_dmuI;
	A(0,3)=de_dmuS;

	A(1,0)=drhoB_dt;
	A(1,1)=drhoB_dmuB;
	A(1,2)=0.5*drhoB_dmuI;
	A(1,3)=drhoB_dmuS;

	A(2,0)=drhoI_dt;
	A(2,1)=drhoI_dmuB;
	A(2,2)=0.5*drhoI_dmuI;
	A(2,3)=drhoI_dmuS;

	A(3,0)=drhoS_dt;
	A(3,1)=drhoS_dmuB;
	A(3,2)=0.5*drhoS_dmuI;
	A(3,3)=drhoS_dmuS;

	//printf("------------------- A(i,j) -------------\n");
	//cout << A << endl;
	//printf("----------------------------------------------\n");
}

double Csampler::GenerateThermalMass(CresInfo *resinfo){
	double mw;
	if(mastersampler->SETMU0)
		mw=maxweight[resinfo->ires];
	else
		mw=maxweight[resinfo->ires];
	double m,m1,m2,kr,k2mr,r1,r2,k,gamma,rho,lor,k2,weight,mass,width;
	bool success;
	double alpha=mastersampler->RESWIDTH_ALPHA;
	CmeanField *mf=mastersampler->meanfield;
	if(resinfo->decay){
		mass=mf->GetMass(resinfo,sigmaf);
		width=resinfo->width;
		m1=resinfo->branchlist[0]->resinfo[0]->mass;
		m1=mf->GetMass(resinfo->branchlist[0]->resinfo[0],sigmaf);
		m2=0.0;
		for(int n=1;n<(resinfo->branchlist[0]->resinfo.size());n++){
			m2+=mf->GetMass(resinfo->branchlist[0]->resinfo[n],sigmaf);
		}
		k2mr = gsl_sf_bessel_Kn(2,(mass/Tf)); // K2 for resmass
		kr=pow(mass*mass-m1*m1-m2*m2,2) - (4*m1*m1*m2*m2);
		kr = (1/(2*mass))*sqrt(abs(kr)); // k at resonant mass
		success=false; // for use in while loop
		do{
			r1 = randy->ran(); // get random numbers
			r2 = randy->ran(); // between [0, 1]
			m = ((width/2)*tan(PI*(r1 - .5))) + mass;// generate random mass value proportional to the lorentz distribution
			if ((m < resinfo->minmass) ) continue;
			// throw out values out of range
        	k=sqrt(abs(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0)))/(2.0*m);
        	if((resinfo->spin)<1.001)
        	{
           	    gamma=width*(mass/m)*(k/kr);
        	}
        	else
        	{
		        gamma=width*(mass/m)*((k*k*k)/(kr*kr*kr))*((kr*kr+HBARC*HBARC)/(k*k+HBARC*HBARC));
        	}
            rho=(2.0)/(width*PI)*0.25*gamma*gamma/((0.25*gamma*gamma)+(mass-m)*(mass-m));

			//k=sqrt(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*m);
			//gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha); // CHANGE??
			//rho=(2.0/(width*PI))*(0.25*gamma*gamma)/((0.25*gamma*gamma)+(mass-m)*(mass-m)); // CHANGE??
            //gamma=width*(resmass/E)*((k*k*k)/(kr*kr*kr))*((kr*kr+HBARC*HBARC)/(k*k+HBARC*HBARC));
            //rho=(2.0*norm/PI)*E*E*gamma/((E*E*gamma*gamma)+(resmass+E)*(resmass+E)*(resmass-E)*(resmass-E));

			lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-m,2.0));
			k2 = gsl_sf_bessel_Kn(2,(m/Tf)); // K2 value
			weight = rho*k2*m*m/(lor*k2mr*mass*mass*mw);
			if (r2 < weight) success=true; // success
		}while(!success);
	}
	else
		m=resinfo->mass;
	return m; //returns a random mass proportional to n0*L'
}

void Csampler::CalcDensitiesF0(){
	CresInfo *resinfo;
	CbranchInfo *bptr;
	CresMassMap::iterator rpos;
	double s,m,m1,m2,degen,width,minmass,code,I3,mutot,mratio=1.0+(sigmaf-93.0)/93.0;
	double Pi,epsiloni,densi,sigma2i,dedti,maxweighti,dm,Tfpi;
	double rhoBcalc=0.0,rhoIcalc=0.0,rhoScalc=0.0;
	int ires,nres,n;
	char dummy[100];
	nhadronsf0=Pf0=epsilonf0=0.0;
	ires=0;
	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
			GetDensPMaxWeight(resinfo,0.0,densi,epsiloni,Pi,dedti,maxweighti);
			densityf0[ires]=densi;
			rhoBcalc+=densi*resinfo->baryon;
			rhoIcalc+=densi*I3;
			rhoScalc+=densi*resinfo->strange;
			maxweight[ires]=maxweighti;
			nhadronsf0+=densityf0[ires];
			epsilonf0i[ires]=epsiloni;
			Pf0i[ires]=Pi;
			epsilonf0+=epsiloni;
			Pf0+=Pi;
			ires+=1;
		}
	}
}

void Csampler::CalcDensitiesF(){
	CresInfo *resinfo;
	CbranchInfo *bptr;
	CresMassMap::iterator rpos;
	double s,m,m1,m2,degen,width,minmass,code,I3,mutot,xx;
	double Pi,epsiloni,densi,sigma2i,dedti,maxweighti,dm;
	double rhoBcalc=0.0,rhoIcalc=0.0,rhoScalc=0.0;
	int ires,nres,n;
	char dummy[100];
	nhadronsf=Pf=epsilonf=0.0;
	ires=0;
	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			GetDensPMaxWeight(resinfo,0.0,densi,epsiloni,Pi,dedti,maxweighti);
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
			mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
			xx=exp(mutot);
			//printf("muB=%g\n",muB);
			//densityf[ires]=densi=densityf0[ires]*xx;
			densityf[ires]=densi=densi*xx;
            //printf("densi=%g\n",densi);
			rhoBcalc+=densi*resinfo->baryon;
			rhoIcalc+=densi*I3;
			rhoScalc+=densi*resinfo->strange;
			nhadronsf+=densi;
			epsilonf+=epsilonf0i[ires]*xx;
			Pf+=Pf0i[ires]*xx;
			ires+=1;
		}
	}
}

void Csampler::GetDensPMaxWeight(CresInfo *resinfo,double mutot,double &densi,double &epsiloni,double &Pi,double &dedti,double &maxweighti){
	double degeni,width,m,minmass,m1,m2,mratio=1.0,sigma2i,dm;
	double ddensi,eepsiloni,PPi,xx;
	int n;
	CmeanField *mf=mastersampler->meanfield;
	//
	bool decay=resinfo->decay;
	//decay=false;
	//
	degeni=resinfo->spin;
	m=mf->GetMass(resinfo,sigmaf);
	if(abs(resinfo->baryon)==1)
		m*=mratio;
	width=resinfo->width;
	minmass=resinfo->minmass;
	if((minmass>0.0) && (width>0.00001) && decay){
        //m1=resinfo->branchlist[0]->resinfo[0]->mass;
		m1=mf->GetMass(resinfo->branchlist[0]->resinfo[0],sigmaf);
		m2=0.0;
		for(n=1;n<(resinfo->branchlist[0]->resinfo.size());n++){
			dm=mf->GetMass(resinfo->branchlist[0]->resinfo[n],sigmaf);
			m2+=dm;
		}
        EOS::freegascalc_onespecies_finitewidth(Tf,m,m1,m2,width,RESWIDTH_ALPHA,degeni,minmass,epsiloni,Pi,densi,sigma2i,dedti,maxweighti);
	}
	else{
		//printf("m=%g, Tf=%g\n",m,Tf);
		EOS::freegascalc_onespecies(Tf,m,epsiloni,Pi,densi,sigma2i,dedti);
		maxweighti=1.0;
        //printf("Tf=%g,m=%g,epsiloni=%g,Pi=%g,densi=%g,sigma2i=%g,dedti=%g\n",Tf,m,epsiloni,Pi,densi,sigma2i,dedti);
	}
    //printf("m=%g, Tf=%g\n",m,Tf);
    //printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
	xx=exp(mutot);
	densi*=degeni*xx;
	Pi*=degeni*xx;
	epsiloni*=degeni*xx;
	dedti*=degeni*xx;
}
