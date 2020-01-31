#ifndef __EOS_H_
#define __EOS_H_
#include "classdefs.h"
#include "resonances.h"

// ------------------------
// functions for calculating EoS (epsilon,P,density,depsilon/dT, and sigma^2) of single species 
// sigma^2 is used for fluctuations 
// freegascalc_onespecies_finitewidth includes spectral information from CresInfo* object
// -----------------------

namespace EOS{
    void freegascalc_onespecies(double T,double m,double &epsilon,double &P,double &dens,double &sigma2,double &dedt);
    void freegascalc_onespecies_finitewidth(CresInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt);
};

#endif
