#ifndef __SAMPLER_DEFS_H__
#define __SAMPLER_DEFS_H__
#include <map>
#include <unordered_map>
#include <vector>
using namespace std;

// -----------------------
// This is a list of class names and typedefs used throughout
// -----------------------

namespace msu_sampler {
  class Cpart;
  class CresList;
  class CresInfo;
  class CbranchInfo;
  class Chyper;
  class CmasterSampler;
  class CmeanField;
  class CdecayInfo;

  typedef map<long int,CresInfo *> CresInfoMap;
  typedef multimap<double,CresInfo *> CresMassMap;
  typedef pair<long int, CresInfo*> CresInfoPair;
  typedef pair<double, CresInfo*> CresMassPair;
  typedef vector<CbranchInfo *> CbranchList; //gives branchlist name
  typedef multimap<long int,Cpart *> CpartMap;
  typedef double FourVector[4];
  typedef map<long int,CdecayInfo *> CdecayInfoMap;
  typedef pair<long int,CdecayInfo *> CdecayInfoPair;
}

#endif
