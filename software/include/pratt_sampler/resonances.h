#ifndef __RESONANCES_H__
#define __RESONANCES_H__
#include <map>
#include <unordered_map>
#include <cmath>
#include <array>
#include <fstream>
#include <gsl/gsl_sf.h>
#include <eigen3/Eigen/Dense>
#include "parametermap.h"
#include "randy.h"
#include "classdefs.h"
#include "misc.h"
#include "constants.h"
using namespace std;

// --------------------
// CresInfo contains all information about an object
// Includes information about decays. Each decay branch is descrbed by CbranchInfo object
// List of resonces is in CresList object
// CresList has 2 maps, resmap has key = pid and value = resinfo ptr, while massmap has key=mass and value = resinfo ptr
// resmap good for finding resinformation indexed by pid, massmap good for thumbing through resonances lightest first
// -------------------

class CbranchInfo{
public:
	vector<CresInfo *> resinfo; //pointers for resinfo
	double branching;
	int L;
	CbranchInfo();
};

class CresInfo{
public:
	int ires;
	double mass;
	double degen;
	double width;
	double minmass;
	string name;
	int pid;
	int charge;
	int strange;
	int baryon;
	int bottom;
	int charm;
	int q[3];
	int up,down;
	int G_Parity;
	bool decay; //false if stable, true if can decay. check if true
	int nchannels;
	CbranchList branchlist;
	CbranchInfo	*bptr_minmass;
	vector<double> SpectVec;
	vector<double> SpectEVec;
	vector<double> GammaVec;
	void Print();
	void PrintBranchInfo();
	void DecayGetResInfoPtr(int &nbodies,array<CresInfo *,5> &daughterresinfo);
	void DecayGetResInfoPtr_minmass(int &nbodies,array<CresInfo *,5> &daughterresinfo);
	bool CheckForDaughters(int pid);
	bool CheckForNeutral();
	double GenerateMass();
	double ChiInt(double T,double vmax); // Integral used by ChiOmega
	double ChiTilde(double T,double vmax); // Integral used by ChiOmega
	CresInfo();
	static Crandy *randy;
	static CresList *reslist;
	static double **ChiA; // stored array used by ChiOmegaInt
	void CalcSpectralFunction();
	void ReadSpectralFunction();
	void CalcMinMass();
	double GetSpectralFunction(double E);
	double GetEofN(int n); // return energy for middle of spectral function bin n
	double GetMeshE(double E); // returns E at middle of mesh cell
	void PrintSpectralFunction();
	double GetRhoAB(double E,CresInfo *resinfo_a,CresInfo *resinfo_b,int L);
	double GetFF(double E,double Ea,double Eb,CresInfo *resinfo_a,CresInfo *resinfo_b);  // form factor for generating spectral functions
	double GetBL2(double k,int L); // factor used to generate spectral functions in arXiv:1606642v2
	double GetBW(double E,double Mr,double Gamma);  // relativistic Breit Wigner
	double GetBW_base(double E,double Mr,double Gamma); // simple non-rel. fixed gamma BW, used as base for Monte Carlo
	double GenerateMass_base();
	double GetDecayMomentum(double M,double ma,double mb);
	void NormalizeSF();  //normalizes spectral function
	bool SFcalculated;
	static string SFDIRNAME; // location of spectral functions for reading in
	static int NSPECTRAL;  // number of points in spectral function
};

class CresList{
public:
	CresList();
	~CresList();
	CresList(CparameterMap* parmap_in);
	CresInfoMap resmap;
	CresMassMap massmap;
	//Cmerge ***MergeArray;
	//double **SigmaMaxArray;
	CresInfo *GetResInfoPtr(int pid);
	void ReadResInfo();
	void CalcMinMasses();
	void ReadSpectralFunctions();
	void CalcSpectralFunctions();
	//bool RESONANCE_DECAYS;
	CparameterMap *parmap;
};

class CdecayInfo{
public:
	int Nparts[25];
	double branchratio[25];
	int products[25][5];
	int d_L[25];
};

#endif
