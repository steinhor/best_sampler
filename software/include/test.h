#ifndef __TEST_H_
#define __TEST_H_

#include "master.h"
using namespace std;

namespace test{
  void n_hyperlist(CmasterSampler ms,int nsamplers);
  void n_dummy(CmasterSampler ms);
  void bose_terms();
  void calc_bose(Csampler *sampler,Chyper *hyper, FILE *file);
  void test_viscous();

};

#endif
