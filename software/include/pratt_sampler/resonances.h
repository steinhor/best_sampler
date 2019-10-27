#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include <map>
#include <unordered_map>
#include <cmath>
#include <array>
#include <fstream>
#include <gsl/gsl_sf.h>
#include <Eigen/Dense>
#include "pratt_sampler/parametermap.h"
#include "pratt_sampler/randy.h"
#include "pratt_sampler/classdefs.h"
#include "pratt_sampler/misc.h"
#include "pratt_sampler/constants.h"

using namespace std;

class CbranchInfo{
public:
	vector<CresInfo *> resinfo; //pointers for resinfo
	double branching;
	CbranchInfo();
};

class CresInfo{
public:
	int ires;
	double mass;
	double spin;
	double width;
	double minmass;
	string name;
	int code;
	int charge;
	int strange;
	int baryon;
	int bottom;
	int charm;
	int q[3];
	int up,down;
	int G_Parity;
	bool decay; //false if stable, true if can decay. check if true
	CbranchList branchlist;
	CbranchInfo	*bptr_minmass;
    map<double,double> spectmap;
	void Print();
	void DecayGetResInfoPtr(int &nbodies,array<CresInfo *,5> &daughterresinfo);
	void DecayGetResInfoPtr_minmass(int &nbodies,array<CresInfo *,5> &daughterresinfo);
	bool CheckForDaughters(int code);
	bool CheckForNeutral();
	double GenerateMass();
	double ChiInt(double T,double vmax); // Integral used by ChiOmega
	double ChiTilde(double T,double vmax); // Integral used by ChiOmega
	CresInfo();
	static Crandy *randy;
	static CresList *reslist;
	static double **ChiA; // stored array used by ChiOmegaInt
};
/*
class Cmerge{
public:
	Cmerge(CresInfo *resinfo,double branching, int L);
	CresInfo *resinfo;
	int L;
	double branching;
	Cmerge *next;
};
*/

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
    void FillSpectMap(CresInfo *resinfo);
	//bool RESONANCE_DECAYS;
	CparameterMap *parmap;
};

#endif
