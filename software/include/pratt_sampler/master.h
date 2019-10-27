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

class CmasterSampler{
public:
	CmasterSampler(CparameterMap *parmapin);
	CmasterSampler(string parsfilename); // Constructor
	~CmasterSampler();
	CparameterMap *parmap;
	CresList *reslist;
	Crandy *randy;
	double TFmin,TFmax;
	double SIGMAFmin,SIGMAFmax;
	int NTF,NSIGMAF;
	double DELTF,DELSIGMAF;
	bool VISCOUSCORRECTIONS,SETMU0;
	double RESWIDTH_ALPHA;
	double YMAX; // only for 2D sampler
	int NEVENTS;
	int nelements;

	list<Chyper *> hyperlist;
	vector<vector<Csampler *>> sampler;  // array of samplers indexed by T and sigma
	vector<Cpart *> part;

	void ReadHyper2D();
	void DeleteVolumeElements();
	int MakeEvent(); // returns number of parts
	Csampler* ChooseSampler(Chyper *hyper);
	void ChooseSampler(double Tf,double sigmaf,int &itf,int &isigmaf);
	void WriteParts();
	static CmeanField *meanfield;
};

#endif
