#ifndef __TEST_H_
#define __TEST_H_

#include "master.h"
using namespace std;

namespace test{
  void n_hyperlist(CmasterSampler ms,int nsamplers);
  void n_dummy(bool readout);
  void bose_terms();
  void calc_bose(Csampler *sampler,Chyper *hyper, FILE *file);
  void test_viscous();
  void newton_method_hyperlist(CmasterSampler ms);
  void newton_method_dummy();
  void eps_rho_derivs();
};

#endif
