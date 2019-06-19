#ifndef __SAMPLER_DEFS_H__
#define __SAMPLER_DEFS_H__

#include <map>
#include <unordered_map>
#include <cmath>
#include <vector>

using namespace std;

class Cpart;
class CresList;
class CresInfo;
class CbranchInfo;
class Chyper;
class CmasterSampler;
class CmeanField;

typedef map<long int,CresInfo *> CresInfoMap;
typedef multimap<double,CresInfo *> CresMassMap;
typedef pair<long int, CresInfo*> CresInfoPair;
typedef pair<double, CresInfo*> CresMassPair;
typedef vector<CbranchInfo *> CbranchList; //gives branchlist name
typedef multimap<long int,Cpart *> CpartMap;
typedef pair<long int,Cpart*> CpartPair;
typedef double FourVector[4];



#endif