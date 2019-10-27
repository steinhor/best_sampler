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
	int n,nbc,imax;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp,J,nfact,sign;
	double lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;

    if (parmap->getB("BOSE_CORR",false)) nbc=parmap->getI("N_BOSE_CORR",1);
    else nbc=1;

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
            if (resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) imax=nbc;
            else imax=1;

            dIpp=0.0;
            for (int i=1;i<=imax;i++) {
    			m=resinfo->mass;
    			degen=resinfo->spin;
    			z=m*i/Tf;
    			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-muS*resinfo->strange);
    			mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
    			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
    			for(int j=1;j<nmax+5;j++){
    				n=5-2*j;
    				if(n!=-1) {
                        G[j]=(-exp(-z)/n)+(G[j-1]*z*z-z*exp(-z))/((n+1.0)*n);
                    }
    				else {
                        G[j]=gsl_sf_gamma_inc(-1,z)*z;
                    }
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
                dIpp+=degen*exp(i*mutot)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
            }
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
	int n,nbc,imax,ires=0;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp,J,nfact,sign,temp;
	double lambdafact,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;

    if (parmap->getB("BOSE_CORR",false)) nbc=parmap->getI("N_BOSE_CORR",1);
    else nbc=1;

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
            if (resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) imax=nbc;
            else imax=1;

            dIpp=0.0;
            for (int i=1;i<=imax;i++) {
    			m=resinfo->mass;
    			degen=resinfo->spin;
    			z=m*i/Tf;

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
                if (resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) {
                    npilambda0[resinfo->code].push_back(temp);
                }
            }
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
	lambdaf=lambdaf0=lambdafact;
}

void Csampler::CalcLambdaF(){
	int n,nbc,ires=0;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp,J,nfact,sign;
	double lambdafact,mutot,I3;
	CresInfo *resinfo;
	CresMassMap::iterator rpos;

    if (parmap->getB("BOSE_CORR",false)) nbc=parmap->getI("N_BOSE_CORR",1);
    else nbc=1;

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-muS*resinfo->strange);
			mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;

            dIpp=0.0;
            if (resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) {
                for (int i=1;i<=nbc;i++) {
                    Ipp+=npilambda0[resinfo->code][i-1]*exp(i*mutot);
                }
            }
            else Ipp+=lambdaf0i[ires]*exp(mutot);
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
    double Pi,epsiloni,densi,dedti,maxweighti;
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

    npiP.clear();
    npiepsilon.clear();
    npidedt.clear();
    npidens.clear();

    for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
        resinfo=rpos->second;
        if(resinfo->code!=22){
            GetDensPMaxWeight(resinfo,densi,epsiloni,Pi,dedti,maxweighti);
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
                if(resinfo->code==211 || resinfo->code==-211) {
                    nh0_b0i2s0+=npidens[resinfo->code][0];
                    eh0_b0i2s0+=npiepsilon[resinfo->code][0];
                    dedth0_b0i2s0+=npidedt[resinfo->code][0];
                }
                else {
                    nh0_b0i2s0+=densi;
                    eh0_b0i2s0+=epsiloni;
                    dedth0_b0i2s0+=dedti;
                }
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

void Csampler::GetMuNH(Chyper *hyper){ //3D Newton's Method
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
    bool bose_corr=parmap->getB("BOSE_CORR",false);
    int n_bose_corr=parmap->getI("N_BOSE_CORR",1);

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
        if (bose_corr) {
            for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
                int n=2;
                if (it->first==211 || it->first==-211) {
                    for(int i=1;i<n_bose_corr;i++) {
                        rhoI+=0.5*((it->second)[i])*(2*exp(muI*n)-2*exp(-muI*n));
                        n++;
                    }
                }
            }
        }

		drhoI_dmuB=drhoB_dmuI;
		drhoI_dmuI=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
			+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);
        if (bose_corr) {
            for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
                int n=2;
                if (it->first==211 || it->first==-211) {
                    for(int i=1;i<n_bose_corr;i++) {
                        drhoI_dmuI+=0.5*((it->second)[i])*(4*exp(muI*n)+4*exp(-muI*n));
                        n++;
                    }
                }
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
		dmu=A.colPivHouseholderQr().solve(drho);

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
    if (bose_corr) {
        for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++) {
                    nhadronsf+=0.5*((it->second)[i])*(exp(muI*n)+exp(-muI*n));
                    n++;
                }
            }
        }
    }
    if (bose_corr) {
        for (CboseMap::iterator it=npiepsilon.begin();it!=npiepsilon.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++) {
                    epsilonf+=0.5*((it->second)[i])*(exp(muI*n)+exp(-muI*n));
                    n++;
                }
            }
        }
    }

	hyper->muB=mu[0];
	hyper->muI=2.0*mu[1];
	hyper->muS=mu[2];
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
	hyper->nhadrons=nhadronsf;
}

void Csampler::GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget){ //4D Newton's Method
	// Here rhoI refers to rho_u-rho_d = 2*I3 and mu[1]=muI/2
    bool bose_corr=parmap->getB("BOSE_CORR",false);
    int n_bose_corr=parmap->getI("N_BOSE_CORR",1);
	double epsilon=0,rhoB=0,rhoS=0,rhoI=0;
	double xB,xI,xS,xxB,xxI,xxS;
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

	}while(fabs(drho[0])>1.0E-3 || fabs(drho[1])>1.0E-5 || fabs(drho[2])>1.0E-5 || fabs(drho[3])>1.0E-5);
	xB=exp(muB);
	xI=exp(0.5*muI);
	xS=exp(muS);
	xxB=1.0/xB;
	xxI=1.0/xI;
	xxS=1.0/xS;
	nhadronsf=nh0_b0i0s0
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
				+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
					+0.5*nh0_b0i2s0*(xI*xI+xxI*xxI)
					+0.5*nh0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
						+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
							+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
								+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
									+0.25*nh0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
    if (bose_corr) {
        for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++) {
                    nhadronsf+=0.5*((it->second)[i])*(exp(muI*n)+exp(-muI*n));
                    n++;
                }
            }
        }
    }
	CalcDensitiesF();
}

void Csampler::GetEpsilonRhoDerivatives(double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A){
	double xB,xI,xS,xxB,xxI,xxS;

	double drhoB_dt,drhoB_dmuB,drhoB_dmuS,drhoB_dmuI;
	double drhoS_dt,drhoS_dmuB,drhoS_dmuS,drhoS_dmuI;
	double drhoI_dt,drhoI_dmuB,drhoI_dmuS,drhoI_dmuI;
	double de_dT,de_dmuB,de_dmuI,de_dmuS;

    bool bose_corr=parmap->getB("BOSE_CORR",false);
    int n_bose_corr=parmap->getI("N_BOSE_CORR",1);

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
    if (bose_corr) {
        for (CboseMap::iterator it=npiepsilon.begin();it!=npiepsilon.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++) {
                    epsilon+=0.5*((it->second)[i])*(exp(muI*n)+exp(-muI*n));
                    n++;
                }
            }
        }
    }

	de_dT=dedth0_b0i0s0+0.5*dedth0_b0i2s0*(xI*xI+xxI*xxI)
		+0.25*dedth0_b0i1s1*(xI+xxI)*(xS+xxS)
			+0.5*dedth0_b1i0s1*(xB*xxS+xxB*xS)
				+0.5*dedth0_b1i0s3*(xB*xxS*xxS*xxS+xxB*xS*xS*xS)
					+0.25*dedth0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*dedth0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*dedth0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI)
								+0.25*dedth0_b1i3s0*(xB+xxB)*(xI*xI*xI+xxI*xxI*xxI);
    if (bose_corr) {
        for (CboseMap::iterator it=npidedt.begin();it!=npidedt.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++){
                    de_dT+=0.5*((it->second)[i])*(exp(muI*n)+exp(-muI*n));
                    n++;
                }
            }
        }
    }

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
    if (bose_corr) {
        for (CboseMap::iterator it=npiepsilon.begin();it!=npiepsilon.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++){
                    de_dmuI+=0.5*((it->second)[i])*(2*exp(muI*n)-2*exp(-muI*n));
                    n++;
                }
            }
        }
    }

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
    if (bose_corr) {
        for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++){
                    rhoI+=0.5*((it->second)[i])*(2*exp(muI*n)-2*exp(-muI*n));
                    n++;
                }
            }
        }
    }

	drhoI_dt=de_dmuI/(Tf*Tf);

	drhoI_dmuB=drhoB_dmuI;

	drhoI_dmuI=0.5*nh0_b0i2s0*(4*xI*xI+4*xxI*xxI)
		+0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
					+0.25*nh0_b1i1s0*(xB+xxB)*(xI+xxI)
						+0.25*nh0_b1i1s2*(xB*xxS*xxS+xxB*xS*xS)*(xI+xxI)
							+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(4*xI*xI+4*xxI*xxI)
								+0.25*nh0_b1i3s0*(xB+xxB)*(9*xI*xI*xI+9*xxI*xxI*xxI);
    if (bose_corr) {
        for (CboseMap::iterator it=npidens.begin();it!=npidens.end();it++) {
            int n=2;
            if (it->first==211 || it->first==-211) {
                for(int i=1;i<n_bose_corr;i++){
                    drhoI_dmuI+=0.5*((it->second)[i])*(4*exp(muI*n)+4*exp(-muI*n));
                    n++;
                }
            }
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

	drhoS_dt=de_dmuS/(Tf*Tf);
	drhoS_dmuB=drhoB_dmuS;
	drhoS_dmuI=drhoI_dmuS;
	drhoS_dmuS=0.25*nh0_b0i1s1*(xI+xxI)*(xS+xxS)
		+0.5*nh0_b1i0s1*(xB*xxS+xxB*xS)
			+0.5*nh0_b1i0s3*(9*xB*xxS*xxS*xxS+9*xxB*xS*xS*xS)
				+0.25*nh0_b1i1s2*(4*xB*xxS*xxS+4*xxB*xS*xS)*(xI+xxI)
					+0.25*nh0_b1i2s1*(xB*xxS+xxB*xS)*(xI*xI+xxI*xxI);

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
}

double Csampler::GenerateThermalMass(CresInfo *resinfo){
    double mw;
    int nfail=0,nfailmax=1000;
    mw=maxweight[resinfo->ires];
    double E,m1,m2,kr,k2mr,r1,r2,k,gamma,rho,lor,k2,weight,mass,width;
    bool success;
    double alpha=mastersampler->RESWIDTH_ALPHA;
    CmeanField *mf=mastersampler->meanfield;
    double decay=resinfo->decay;
    decay=false;
    if(decay){
        mass=mf->GetMass(resinfo,sigmaf);
        width=resinfo->width;
        k2mr = gsl_sf_bessel_Kn(2,(mass/Tf)); // K2 for resmass
        success=true; // for use in while loop
        do{
            r1 = randy->ran(); // get random numbers
            r2 = randy->ran(); // between [0, 1]
            E = ((width/2)*tan(PI*(r1 - .5))) + mass;// generate random mass value proportional to the lorentz distribution
            if ((E < resinfo->minmass) ) continue;
            k2 = gsl_sf_bessel_Kn(2,(E/Tf)); // K2 value
            auto it = resinfo->spectmap.lower_bound(E);
            rho = (*it).second;
            lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-E,2.0));
            rho = lor;
            weight = rho*k2*E*E/(lor*k2mr*mass*mass*mw);
            if (r2 < weight) success=true; // success
            else nfail+=1;
        }while(!success && nfail<nfailmax);
    }
    else
        E=resinfo->mass;
    return E; //returns a random mass proportional to n0*L'
}

void Csampler::CalcDensitiesF0(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double I3,Pi,epsiloni,densi,dedti,maxweighti;
	int ires=0;
	nhadronsf0=Pf0=epsilonf0=0.0;

    npiP.clear();
    npiepsilon.clear();
    npidedt.clear();
    npidens.clear();

    npiP0.clear();
    npiepsilon0.clear();
    npidens0.clear();

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
			GetDensPMaxWeight(resinfo,densi,epsiloni,Pi,dedti,maxweighti);
			densityf0[ires]=densi;
            epsilonf0i[ires]=epsiloni;
			Pf0i[ires]=Pi;

            if ((resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) && parmap->getB("BOSE_CORR",false)) {
                for (int i=0;i<parmap->getI("N_BOSE_CORR",1);i++) {
                    npidens0[resinfo->code].push_back(npidens[resinfo->code][i]);
                    npiepsilon0[resinfo->code].push_back(npiepsilon[resinfo->code][i]);
                    npiP0[resinfo->code].push_back(npiP[resinfo->code][i]);
                }
            }

			maxweight[ires]=maxweighti;
			nhadronsf0+=densityf0[ires];
			epsilonf0+=epsiloni;
			Pf0+=Pi;
			ires+=1;
		}
	}
}

void Csampler::CalcDensitiesF(){
	CresInfo *resinfo;
	CresMassMap::iterator rpos;
	double I3,mutot,xx,Pi,epsiloni,densi,dedti,maxweighti,dm;
    int ires=0;
	nhadronsf=Pf=epsilonf=0.0;

    npiP.clear();
    npiepsilon.clear();
    npidedt.clear();
    npidens.clear();

	for(rpos=reslist->massmap.begin();rpos!=reslist->massmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
            I3=0.5*(2.0*resinfo->charge-resinfo->baryon-resinfo->strange);
            mutot=muB*resinfo->baryon+muI*I3+muS*resinfo->strange;
			xx=exp(mutot);
            GetDensPMaxWeight(resinfo,densi,epsiloni,Pi,dedti,maxweighti);
            if((resinfo->code==211 || resinfo->code==-211 || resinfo->code==111) && parmap->getB("BOSE_CORR",false)) {
                densi=0;
                for (int i=0;i<parmap->getI("N_BOSE_CORR",1);i++) {
                    densi+=npidens[resinfo->code][i]*exp(mutot*(i+1));
                    epsilonf+=npiepsilon[resinfo->code][i]*exp(mutot*(i+1));
                    Pf+=npiP[resinfo->code][i]*exp(mutot*(i+1));
                }
            }
            else {
                densi=densi*xx;
                epsilonf+=epsiloni*xx;
                Pf+=Pi*xx;
            }
            densityf[ires]=densi;
            nhadronsf+=densi;
            ires++;
		}
	}
}

void Csampler::GetDensPMaxWeight(CresInfo *resinfo,double &densi,double &epsiloni,double &Pi,double &dedti,double &maxweighti){
	double degen,width,m,minmass,m1,m2,sigma2i,dm;
	CmeanField *mf=mastersampler->meanfield;
	bool decay=resinfo->decay;
	//decay=false;
	degen=resinfo->spin;
	m=mf->GetMass(resinfo,sigmaf);
	width=resinfo->width;
	minmass=resinfo->minmass;

  if(Tf==0) {
    printf("%d",resinfo->code);
    exit(1);
  }

	if((minmass>0.0) && (width>0.00001) && decay){
		m1=mf->GetMass(resinfo->branchlist[0]->resinfo[0],sigmaf);
		m2=0.0;
		for(int n=1;n<(resinfo->branchlist[0]->resinfo.size());n++){
			dm=mf->GetMass(resinfo->branchlist[0]->resinfo[n],sigmaf);
			m2+=dm;
		}
        EOS::freegascalc_onespecies_finitewidth(npidens,npiP,npiepsilon,npidedt,parmap,resinfo,Tf,m,m1,m2,RESWIDTH_ALPHA,epsiloni,Pi,densi,sigma2i,dedti,maxweighti);
    }
	else{
		EOS::freegascalc_onespecies(npidens,npiP,npiepsilon,npidedt,parmap,resinfo,Tf,m,epsiloni,Pi,densi,sigma2i,dedti);
		maxweighti=1.0;
	}
}
