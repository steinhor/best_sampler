#ifndef __SAMPLER_MEANFIELD_H__
#define __SAMPLER_MEANFIELD_H__

#include "classdefs.h"
#include "resonances.h"
#include "parametermap.h"

namespace msu_sampler {
	class CmeanField{
	public:
		bool crap;
		CmeanField();
		virtual double GetMass(CresInfo* resinfo,double sigma);
	};

	class CmeanField_Simple : public CmeanField {
	public:
		CmeanField_Simple(CparameterMap *parmap);
		double GetMass(CresInfo* resinfo,double sigma);
	};
}

#endif
