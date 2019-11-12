#ifndef __TEST_H_
#define __TEST_H_

#include "master.h"
using namespace std;

namespace test{
  void n_hyperlist(CmasterSampler ms,int nsamplers);
  void n_dummy(bool readout);
  void bose_terms(CmasterSampler *ms);
  void calc_bose(CmasterSampler *ms,Csampler *sampler,Chyper *hyper, FILE *file);
  void newton_method_hyperlist(CmasterSampler ms);
  void newton_method_dummy();
};

#endif
