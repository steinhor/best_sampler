#ifndef __SAMPLER_H__
#define __SAMPLER_H__
#include "resonances.h"
#include "part.h"
#include <Eigen/Dense>
#include "classdefs.h"
#include "master.h"
using namespace std;

// -----------------
// This is the workhorse class
// This stores information and methods required to sample particles, but not the specific hyper-element information
// Each sampler object has a specific temperature and sigma field (sigma field considered constant currently)
// If two hyper-elements have same T but different mu, sampler adjusts on the fly to account for new mu
// but stored coefficients that depend on T are common
// -----------------

namespace msu_sampler {
	class Csampler{
	public:
		// Vary for each Csampler instance
		double Tf,sigmaf;
		// functions of Tf,sigma, mu
		// Same but with mu=0
		double nhadrons0,epsilon0,P0,p4overE30;
		vector<double> density0i,epsilon0i,P0i;
		vector<map<double,double>> sfdens0imap;
		bool FIRSTCALL;
		int ntest;

		double GenerateThermalMass(CresInfo *resinfo);
		void GetDensPMax(CresInfo *resinfo,double &densi,double &epsiloni,double &Pi,double &dedti,double &p4overE3);

		Csampler(double Tfset,double sigmafset); // Constructor
		~Csampler();
		void GetPars(CparameterMap *parmapset);
		CpartMap *partmap;
		void CalcLambdaMu0(); // Calculates lambda which is used for viscous corrections with mu=0
		double CalcLambdaF(Chyper *hyper);
		double CalcLambdaF(double muB,double muI,double muS,double Pf); //  calculates lambda with mu !=0
		void CalcDensitiesF();
		void CalcDensitiesMu0();

		double totvol;
		vector<double> pibose_dens0,pibose_P0,pibose_epsilon0,pibose_dedt0;

		//some variables used only for testing
		vector<double> dN_dp_p2;
		double dp;
		bool bose_test_off;
		bool bose_test;
		bool viscous_test;
		bool byflavor_calculated;
		map<int,vector<Cpart*>> pmap;
		map<int,double> DensityMap;

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
		// p^4/E^3, used for viscosity 
		double p4overE3h0_b0i0s0,p4overE3h0_b0i2s0,p4overE3h0_b0i1s1;
		double p4overE3h0_b1i0s1,p4overE3h0_b1i0s3;
		double p4overE3h0_b1i1s0,p4overE3h0_b1i1s2;
		double p4overE3h0_b1i2s1,p4overE3h0_b1i3s0;
		double p4overE3h0_b2i0s0;
		// epslion_h^2/n_h
		double eEbarh0_b0i0s0,eEbarh0_b0i2s0,eEbarh0_b0i1s1;
		double eEbarh0_b1i0s1,eEbarh0_b1i0s3;
		double eEbarh0_b1i1s0,eEbarh0_b1i1s2;
		double eEbarh0_b1i2s1,eEbarh0_b1i3s0;
		double eEbarh0_b2i0s0;
		// m^2*dens
		double m2densh0_b0i0s0,m2densh0_b0i2s0,m2densh0_b0i1s1;
		double m2densh0_b1i0s1,m2densh0_b1i0s3;
		double m2densh0_b1i1s0,m2densh0_b1i1s2;
		double m2densh0_b1i2s1,m2densh0_b1i3s0;
		double m2densh0_b2i0s0;

		void GetNHMu0(); // Calculates above quantities
		void GetMuNH(Chyper *hyper);
		void GetMuNH(double rhoBtarget,double rhoItarget,double rhoStarget,double &muB,double &muI,double &muS);
		void CalcNHadronsEpsilonP(Chyper *hyper);
		void CalcNHadronsEpsilonP(double muB,double muI,double muS,double &nhadronsf,double &epsilonf,double &Pf);
		void GetTfMuNH(Chyper *hyper);
		void GetTfMuNH(double epsilontarget,double rhoBtarget,double rhoItarget,double rhoStarget,double &muB,double &muI,double &muS);
		void GetEpsilonRhoDerivatives(double muB,double muI,double muS,double &epsilon,double &rhoB,double &rhoI,double &rhoS,Eigen::MatrixXd &A);
		int MakeParts(Chyper *hyper);
		void CalcRvisc(Chyper *hyper);
		void BulkScale(Chyper *hyper,double mass,FourVector &pnobulk,FourVector &p);
		void ShearScale(Chyper *hyper,double mass,FourVector &pnoshear,FourVector &p);
		int CheckResInVolume(double dN,double T,CresInfo *resinfo,Chyper *hyper);
		void GetP(Chyper *hyper,double T,CresInfo *resinfo,FourVector &p);
		void GetPInFluidFrame(double m,Chyper *hyper,double T,FourVector &p);
		void CalcSFDensMap(CresInfo *resinfo,double T,map<double,double> &sfdensmap);

		static Crandy *randy;
		static CresList *reslist;
		static CparameterMap *parmap;
		static CmasterSampler *mastersampler;
		static bool CALCMU;
		static bool bose_corr;
		static int n_bose_corr;
		static bool BJORKEN_2D;  // these variables are only used for Bjorken
		static double BJORKEN_YMAX;
		static bool USE_POLE_MASS;
		static bool INCLUDE_BULK_VISCOSITY;
		static bool INCLUDE_SHEAR_VISCOSITY;
	};
}

#endif
