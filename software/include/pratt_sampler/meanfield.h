#ifndef __SAMPLER_MEANFIELD_H__
#define __SAMPLER_MEANFIELD_H__

#include <cmath>
#include "classdefs.h"
#include "resonances.h"
#include "parametermap.h"

namespace pratt_sampler {
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
