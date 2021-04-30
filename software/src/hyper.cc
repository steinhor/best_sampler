#include "msu_sampler/hyper.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "msu_sampler/sampler.h"

using namespace msu_sampler;

Chyper::Chyper(){
	muB=muI=muS=0.0;
	Rvisc_calculated=false;
	firstcall=true;
	int alpha,beta;
	sigma=0.093;
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++)
			pitilde[alpha][beta]=0.0;
	}
}

void Chyper::SetSampler(Csampler *samplerptr){
	sampler=samplerptr;
}

void Chyper::Print(){
	printf("----- hyper info -------\n");
	printf("T0=%g, sampler->Tf=%g, sigma=%g, rhoB=%g, rhoI=%g\n",T0,sampler->Tf,sigma,rhoB,rhoI);
	printf("muB=%g, muI=%g, muS=%g\n",muB,muI,muS);
	printf("epsilon=%g, P=%g, s=%g\n",epsilon,P,GetEntropyDensity());
	printf("dOmega=(%g,%g,%g,%g), udotdOmega=%g\n",dOmega[0],dOmega[1],dOmega[2],dOmega[3],udotdOmega);
	printf("u=(%g,%g,%g,%g).    tau=%g,eta=%g\n",u[0],u[1],u[2],u[3],tau,eta);
	for(int alpha=0;alpha<4;alpha++){
			printf("%10.3e %10.3e %10.3e %10.3e\n",
			pitilde[alpha][0],pitilde[alpha][1],pitilde[alpha][2],pitilde[alpha][3]);
	}
	printf("-----------------------\n");
}

double Chyper::GetEntropyDensity(){
	double s=((epsilon+P)/T0)-muB*rhoB-muI*rhoI-muS*rhoS;
	return s;
}

void Chyper::CalcBiggestpitilde(){
	int alpha,beta;
	Eigen::MatrixXd A(3,3);
	Eigen::VectorXd V(3);
	A(0,0)=0.0;
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++)
			A(alpha-1,beta-1)=pitilde[alpha][beta];
	}
	//Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
	cout << A << endl;
	cout << es.eigenvalues() << endl;
	V=es.eigenvalues();
	printf("eigenvalues: %g,%g,%g\n",V(0),V(1),V(2));
	biggestpitilde=0.0;
	for(alpha=0;alpha<3;alpha++){
		if(V(alpha)>biggestpitilde)
			biggestpitilde=V(alpha);
	}
}
