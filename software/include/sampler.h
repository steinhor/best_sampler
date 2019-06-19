#ifndef __SAMPLER_H__
#define __SAMPLER_H__
#include "resonances.h"
#include "part.h"
#include <eigen3/Eigen/Dense>
#include "classdefs.h"
#include "master.h"

using namespace std;

class Csampler{
public:
	// Vary for each Csampler instance
	double Tf,sigmaf;
	// functions of Tf,sigma, mu
	double muB,muI,muS;
	double nhadronsf,epsilonf,Pf;
	vector<double> densityf;
	vector<double>piplusdensf,piminusdensf,pi0densf;
	double lambdaf; // Functions of T & sigma used for viscous corrections
	// Same but with mu=0
	double nhadronsf0,epsilonf0,Pf0,lambdaf0;
	vector<double> densityf0,epsilonf0i,Pf0i,lambdaf0i;
	vector<double> pidensf0,pi0densf0;
	//
	vector<double> maxweight;
	double GenerateThermalMass(CresInfo *resinfo);
	void GetDensPMaxWeight(CresInfo *resinfo,double mutot,double &densi,double &epsiloni,double &Pi,double &dedti,double &maxweighti);

	bool VISCOUSCORRECTIONS;
	Csampler(double Tfset,double sigmafset); // Constructor
	~Csampler();
	void GetPars(CparameterMap *parmapset);
	CpartMap *partmap;
	double ETAMAX; // |maximum |spatial rapidity|
	void CalcLambda(); // Deprecated
	void CalcLambdaF0(); // Calculates lambda which is used for viscous corrections with mu=0
	void CalcLambdaF(); //  calculates lambda with mu !=0

	void CalcDensitiesF();

	void CalcDensitiesF0();

	// Including Isospin (i refers to 2*I3)
	// numberdensities
	double nh0_b0i0s0,nh0_b0i2s0,nh0_b0i1s1;
	double nh0_b1i0s1,nh0_b1i0s3;
	double nh0_b1i1s0,nh0_b1i1s2;
	double nh0_b1i2s1,nh0_b1i3s0;
    double nh0_b2i0s0;
	// energy densities
	double eh0_b0i0s0,eh0_b0i2s0,eh0_b0i1s1;
	double eh0_b1i0s1,eh0_b1i0s3;
	double eh0_b1i1s0,eh0_b1i1s2;
	double eh0_b1i2s1,eh0_b1i3s0;
    double eh0_b2i0s0;
	// depsilon/dT
	double dedth0_b0i0s0,dedth0_b0i2s0,dedth0_b0i1s1;
	double dedth0_b1i0s1,dedth0_b1i0s3;
	double dedth0_b1i1s0,dedth0_b1i1s2;
	double dedth0_b1i2s1,dedth0_b1i3s0;
    double dedth0_b2i0s0;

	void GetNH0(); // Calculates above quantities
	void GetMuNH(Chyper *hyper);
	void GetTfMuNH(Chyper *hyper);
	void GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget);
	void GetEpsilonRhoDerivatives(double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A);
	int MakeParts(Chyper *hyper);
	void GetP(Chyper *hyper,CresInfo *resinfo,FourVector &p);

	// Same as above with with rhoI != 0, I1 refers to |I_3|=1/2, I3 to 3/2...
	void GetMuNH(double rhoBtarget,double rhoItarget,double &muB,double &muI,double &muS,double &nh); // Uses nh0_xxx, rhoB, rhoI and rhoS=0 to find muB/T, muI/T and muS/T, also finds total hadron density
	double RESWIDTH_ALPHA;
	bool RESONANCE_DECAYS;
	bool USEPOLEMASS;

	static Crandy *randy;
	static CresList *reslist;
	static CparameterMap *parmap;
	static CmasterSampler *mastersampler;
};



#endif
