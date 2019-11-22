#ifndef __MASTER_SAMPLER_H__
#define __MASTER_SAMPLER_H__
#include "classdefs.h"
#include "meanfield.h"
#include "constants.h"
#include "parametermap.h"
#include "randy.h"
#include "misc.h"
#include "resonances.h"
#include <list>
#include <iterator>
//#include "sampler.h"
#include "hyper.h"
#include "sampler.h"
#include "eos.h"
using namespace std;

// --------------------------------
// This object contains array of sampler objects, and information used across all samplers, e.g. resonance info lists
// sampler objects have unique temperature and sigma field
// It also has list of hyper-volume elements
//
// Typical Usage:
// CmasterSampler ms(parametermap);
// ms.partlist=new CpartList(&parmap);
// ms.ReadHyper2D();
// nparts+=ms.MakeEvent();
//
// ---------------------------------

class CmasterSampler{
public:
	CmasterSampler(CparameterMap *parmapin);
	CmasterSampler(string parsfilename); // Constructor
	~CmasterSampler();
	CparameterMap *parmap;
	CresList *reslist;
	Crandy *randy;
	double TFmin,TFmax;
	string MEANFIELD;
	double SIGMAFmin,SIGMAFmax;
	int NTF,NSIGMAF;
	double DELTF,DELSIGMAF;
	bool VISCOUSCORRECTIONS,SETMU0,CALCMU;
	double RESWIDTH_ALPHA;
	double YMAX; // only for 2D sampler
	int NEVENTS,NEVENTS_TOT;
	int nelements;
	bool FINDT;  // true if you need to find T(epsilon) in hyper elements

	list<Chyper *> hyperlist;
	vector<vector<Csampler *>> sampler;  // array of samplers indexed by T and sigma
	CpartList *partlist;

	void ReadHyper2D();
	int MakeEvent(); // returns number of parts
	Csampler* ChooseSampler(Chyper *hyper);
	void ChooseSampler(double Tf,double sigmaf,int &itf,int &isigmaf);
	void MakeDummyHyper();
	static CmeanField *meanfield;
};


#endif
