#ifndef __PART_H__
#define __PART_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <map>
#include <unordered_map>
#include <sys/stat.h>
#include "misc.h"
#include "eos.h"
#include "resonances.h"
#include "classdefs.h"

using namespace std;
//using namespace boost;

class Cpart{
public:
	Cpart();
	~Cpart();
	double tau0;
	double y,eta,msquared;
	FourVector p,r;
 	CresInfo *resinfo;
	
	void Init(int pid,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity);
	void Setp0();
	void SetMass();
	void Copy(Cpart *part);
	void CopyPositionInfo(Cpart *part);
	void CopyMomentumInfo(Cpart *part);
	double GetMass();
	void SetInitialKey();
	void Kill();
	void Print();
	void BjorkenTranslate();
	void BjorkenUnTranslate();
	void Boost(FourVector &u);
	void BoostP(FourVector &u);
	void BoostR(FourVector &u);
	double GetEta(double tau);
	double GetPseudoRapidity();
	double GetRapidity();
	double GetMT();
	void SetY();
	void SetEta(double neweta);
	void GetHBTPars(double &t,double &rout,double &rside,double &rlong);
	CpartMap::iterator GetPos(CpartMap *pmap);
	CpartMap::iterator DeleteFromMap(CpartMap *partmap);
	static CresList *reslist;
};

#endif
