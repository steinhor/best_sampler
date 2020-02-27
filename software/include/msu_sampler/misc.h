#ifndef __INCLUDE_SAMPLER_MISC_H_
#define __INCLUDE_SAMPLER_MISC_H_
#define __MISC_USE_GSL__ // use gnu sci. lib. in misc.cc
#include "classdefs.h"
using namespace std;

// --------------------------------
// Some random functions -- not all are used, or are copied into the BEST directory
//

namespace msu_sampler {
	namespace Misc{
		void BoostToCM(FourVector &u,FourVector &p,FourVector &ptilde);
		void Boost(FourVector &u,FourVector &ptilde,FourVector &p);

		void Pause();
		void Pause(int seconds);
	};
}

#endif
