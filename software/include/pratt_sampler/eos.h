#ifndef __EOS_H_
#define __EOS_H_
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <boost/math/special_functions/bessel.hpp>
#include "classdefs.h"
#include "constants.h"
#include "sf.h"
#include "resonances.h"
using namespace std;

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
