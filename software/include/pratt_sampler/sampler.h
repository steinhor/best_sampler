#ifndef __SAMPLER_H__
#define __SAMPLER_H__
#include "pratt_sampler/resonances.h"
#include "pratt_sampler/part.h"
#include <Eigen/Dense>
#include "pratt_sampler/classdefs.h"
#include "pratt_sampler/master.h"

using namespace std;

class Csampler{
public:
	// Vary for each Csampler instance
	double Tf,sigmaf;
	// functions of Tf,sigma, mu
	double muB,muI,muS;
	double nhadronsf,epsilonf,Pf;
	vector<double> densityf;
	// Same but with mu=0
	double nhadronsf0,epsilonf0,Pf0,lambdaf0,lambdaf;
	vector<double> densityf0,epsilonf0i,Pf0i,lambdaf0i;

	vector<double> maxweight;
	double GenerateThermalMass(CresInfo *resinfo);
	void GetDensPMaxWeight(CresInfo *resinfo,double &densi,double &epsiloni,double &Pi,double &dedti,double &maxweighti);

	bool VISCOUSCORRECTIONS;
	Csampler(double Tfset,double sigmafset); // Constructor
	~Csampler();
	void GetPars(CparameterMap *parmapset);
	CpartMap *partmap;
	void CalcLambda(); // Deprecated
	void CalcLambdaF0(); // Calculates lambda which is used for viscous corrections with mu=0
	void CalcLambdaF(); //  calculates lambda with mu !=0

	void CalcDensitiesF();
	void CalcDensitiesF0();

	map<int,double> DensityMap;
	double totvol;

	//some variables used only for testing
	vector<double> dN_dp_p2;
	double dp;
	bool bose_test_off;
	bool bose_test;
	bool viscous_test;
	map<int,vector<Cpart*>> pmap;

	CboseMap npidens, npiP, npiepsilon, npidedt;
	CboseMap npidens0, npiP0, npiepsilon0, npilambda0;
	bool bose_corr;
	int n_bose_corr;

	// Including Isospin (i refers to 2*I3)
	// number densities
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
	void GetNHadronsf(Chyper *hyper);
	void GetTfMuNH(Chyper *hyper);
	void GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget);
	void GetEpsilonRhoDerivatives(double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A);
	int MakeParts(Chyper *hyper);
	void GetP(Chyper *hyper,CresInfo *resinfo,Cpart *part,double T);

	// Same as above with with rhoI != 0, I1 refers to |I_3|=1/2, I3 to 3/2...
	void GetMuNH(double rhoBtarget,double rhoItarget,double &muB,double &muI,double &muS,double &nh); // Uses nh0_xxx, rhoB, rhoI and rhoS=0 to find muB/T, muI/T and muS/T, also finds total hadron density
	double RESWIDTH_ALPHA;

	static Crandy *randy;
	static CresList *reslist;
	static CparameterMap *parmap;
	static CmasterSampler *mastersampler;
	static bool CALCMU;
};

#endif
