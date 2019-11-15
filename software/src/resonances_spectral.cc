#ifndef __RESONANCES_SPECTRAL_CC__
#define __RESONANCES_SPECTRAL_CC__
#include "pratt_sampler/resonances.h"

void CresList::CalcSpectralFunctions(){
	CresMassMap::iterator rpos;
	CresInfo *resinfo;
	for(rpos=massmap.begin();rpos!=massmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->decay && !resinfo->SFcalculated){
			resinfo->CalcSpectralFunction();
		}
		else
			resinfo->SFcalculated=true;
	}
}

void CresInfo::CalcSpectralFunction(){
	if(!decay){
		printf("Calling CalcSpectralFunction() for stable particle\n");
		exit(1);
	}
	int n,ibranch;
	double E,M0=mass,A,rhoratio;
	double rho_ab,rho_ab0,Gamma,Gamma0;
	CresInfo *resinfo_a,*resinfo_b;
	CresInfo *resinfoswitch;
	if(spectvec.size()==0)
		spectvec.resize(NSPECTRAL);
	for(n=0;n<NSPECTRAL;n++){
		E=GetEofN(n);
		Gamma=0;
		//if(pid==229 && E>minmass)
			//printf("E=%6.4f ",E);
		for(ibranch=0;ibranch<branchlist.size();ibranch++){
			Gamma0=width*branchlist[ibranch]->branching;
			resinfo_a=branchlist[ibranch]->resinfo[0];
			resinfo_b=branchlist[ibranch]->resinfo[1];
			if(E>resinfo_a->minmass+resinfo_b->minmass){
				if(resinfo_b->decay && !resinfo_a->decay){
					resinfoswitch=resinfo_a;
					resinfo_a=resinfo_b;
					resinfo_b=resinfoswitch;
				}
				rho_ab0=GetRhoAB(mass,resinfo_a,resinfo_b,branchlist[ibranch]->L);
				if(E>minmass && rho_ab0>1.0E-5){
					rho_ab=GetRhoAB(E,resinfo_a,resinfo_b,branchlist[ibranch]->L);
					rhoratio=rho_ab/rho_ab0;
					//if(pid==229)
						//printf("%8.4f ",rhoratio);
					/*
					if(rhoratio>300.0){
						printf("--------------------------------------\n");
						printf("spectral function correction = %g\n",rhoratio);
						printf("pid=%d, mass=%g, E=%g, width=%g, branching=%g\n",pid,mass,E,width,branchlist[ibranch]->branching);
						printf("masses of products=%g,%g\n",resinfo_a->mass,resinfo_b->mass);
						printf("pf=%g, pf0=%g, (E-M0)/Gamma=%g\n",GetDecayMomentum(E,resinfo_a->mass,resinfo_b->mass),GetDecayMomentum(mass,resinfo_a->mass,resinfo_b->mass),(E-mass)/width);
					}*/
					Gamma+=Gamma0*rhoratio;
					if(Gamma!=Gamma){
						printf("Disaster, Gamma=%g\n",Gamma);
						Print();
						printf("---------------------------------\n");
						printf("ibranch=%d, E=%g, Gamma_ab=%g\n",ibranch,E,Gamma0*rho_ab/rho_ab0);
						printf("Gamma0=%g, rho_ab=%g, rho_ab0=%g\n",Gamma0,rho_ab,rho_ab0);
						printf("mass-m_amin-m_bmin=%g\n",mass-resinfo_a->minmass-resinfo_b->minmass);
						resinfo_a->Print();
						resinfo_b->Print();
						exit(1);
					}
				}
			}
		}
		//if(pid==229 && E>minmass)
		//	printf("\n");
		A=GetBW(E,M0,Gamma);
		spectvec[n]=A/GetBW_base(E,M0,Gamma0);
	}
	NormalizeSF();
	SFcalculated=true;
}

double CresInfo::GetSpectralFunction(double E){
	int n;
	n=floorl(NSPECTRAL*(0.5+atan2(E-mass,0.5*width)/PI));
	return spectvec[n];
}

double CresInfo::GetEofN(int n){
	return mass+0.5*width*tan(PI*(double(n)+0.5-0.5*NSPECTRAL)/double(NSPECTRAL));
}

double CresInfo::GetMeshE(double E){
	int n=floorl(NSPECTRAL*(0.5+atan2(E-mass,0.5*width)/PI));
	return GetEofN(n);
}

double CresInfo::GetRhoAB(double E,CresInfo *resinfo_a,CresInfo *resinfo_b,int L){
	//printf("CHECK, starting CresInfo::GetRhoAB()\n");
	double pf=0.0,rho_ab=0.0,drho=0.0;
	double Ea,Aa,Eb,Ab;
	int na,nb;
	map<double,double>::iterator ita,itb;
	if(!resinfo_a->decay){  // both a & b stable
		pf=GetDecayMomentum(E,resinfo_a->mass,resinfo_b->mass);
		rho_ab=(pf/E)*GetBL2(pf,L)*GetFF(E,resinfo_a->mass,resinfo_b->mass,resinfo_a,resinfo_b);
	}
	else if(!resinfo_b->decay){ // a is unstable, b is stable
		if(!(resinfo_a->SFcalculated))
			resinfo_a->CalcSpectralFunction();
		for(na=0;na<NSPECTRAL;na++){
			Ea=resinfo_a->GetEofN(na);
			if(Ea>=resinfo_a->minmass){
				Aa=resinfo_a->spectvec[na];
				pf=GetDecayMomentum(E,Ea,resinfo_b->mass);
				drho=Aa*(pf/E)*GetBL2(pf,L)*GetFF(E,Ea,resinfo_b->mass,resinfo_a,resinfo_b);
				rho_ab+=drho;
			}
		}
	}
	else{   // both a and b are unstable
		if(!(resinfo_a->SFcalculated))
			resinfo_a->CalcSpectralFunction();
		if(!(resinfo_b->SFcalculated))
			resinfo_b->CalcSpectralFunction();
		for(na=0;na<NSPECTRAL;na++){
			Ea=resinfo_a->GetEofN(na);
			if(Ea>=resinfo_a->minmass){
				Aa=resinfo_a->spectvec[na];
				for(nb=0;nb<NSPECTRAL;nb++){
					Eb=resinfo_b->GetEofN(nb);
					if(Ea+Eb<E && Eb>=resinfo_b->minmass){
						Ab=resinfo_b->spectvec[nb];
						pf=GetDecayMomentum(E,Ea,Eb);
						drho=Aa*Ab*(pf/E)*GetBL2(pf,L)*GetFF(E,Ea,Eb,resinfo_a,resinfo_b);
						rho_ab+=drho;
					}
				}
			}
		}
	}
	return rho_ab;
}

double CresInfo::GetFF(double E,double Ea,double Eb,CresInfo *resinfo_a,CresInfo *resinfo_b){
	// For Factor From Eq. (36), arXiv:1606.06642v2
	double lambda,s0,FF;
	s0=(Ea+Eb)*(Ea+Eb);   // really? not on-shell masses?

	if(resinfo_a->decay && resinfo_b->decay){
		lambda=0.6;
	}
	else if(baryon!=0 && (resinfo_a->decay || resinfo_b->decay)){
		lambda=2.0;
	}
	else if(baryon==0 && (resinfo_a->decay || resinfo_b->decay)){
		lambda=1.6;
	}
	if(pid==223)
		lambda=0.8;

	FF=(pow(lambda,4)+0.25*pow(s0-mass*mass,2))
		/(pow(lambda,4)+pow( E*E-0.5*(s0+mass*mass),2));
	FF=1.0;
	return FF*FF;
}

double CresInfo::GetBL2(double k,int L){
	const double R=1.0;
	double BL2;
	double x2;
	if(L==0)
		return 1.0;
	else{
		x2=pow(k*R/HBARC,2);
		BL2=pow(x2/(1.0+x2),L);
	}
	return BL2;
}

double CresInfo::GetBW(double E,double M0,double Gamma){
	return (2.0*E*E*Gamma/PI)/(pow(E*E-M0*M0,2)+E*E*Gamma*Gamma);
}

double CresInfo::GetBW_base(double E,double M0,double Gamma0){ // simple lorentzian used for MC weights
	return (0.5*Gamma0/PI)/((E-M0)*(E-M0)+0.25*Gamma0*Gamma0);
	//return (2.0*E*E*Gamma0/PI)/(pow(E*E-M0*M0,2)+E*E*Gamma0*Gamma0);
}

double CresInfo::GenerateMass_base(){
	double r1=randy->ran();
	return ((width/2)*tan(PI*(r1-0.5))) + mass;  // mass according to BW distribution
}

void CresInfo::NormalizeSF(){
	int n;
	double norm=0.0;
	for(n=0;n<NSPECTRAL;n++){
		norm+=spectvec[n];
	}
	for(n=0;n<NSPECTRAL;n++){
		spectvec[n]*=NSPECTRAL/norm;
	}
}

void CresInfo::PrintSpectralFunction(){
	int n;
	double E,k2,k2mr,Tf=0.150,basedist;
	printf("- Spectral Function for pid=%d -\n",pid);
	printf("__ E ___ A/A_BW ___ A _ (A/A_BW)*dens(M)/dens(M0) _ A*dens(M)/dens(M0) \n");
	k2mr=mass*mass*Bessel::Kn(2,mass/Tf);
	for(n=0;n<NSPECTRAL;n++){
		E=GetEofN(n);
		k2=E*E*Bessel::Kn(2,E/Tf);
		basedist=GetBW_base(E,mass,width);
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f\n",
		E,spectvec[n],spectvec[n]*basedist,spectvec[n]*k2/k2mr,spectvec[n]*k2*basedist/k2mr);
	}
}


double CresInfo::GetDecayMomentum(double M,double ma,double mb){ // Gives relative momentum
	if(ma+mb<M)
		return sqrt(fabs(pow((M*M-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*M);
	else
		return 0.0;
}

#endif
