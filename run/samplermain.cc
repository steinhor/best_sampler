#include "msu_sampler/master.h"
#include <cstring>
using namespace std;
using namespace msu_sampler;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	CpartList pl=CpartList(&parmap, ms.reslist);
	ms.partlist=&pl;
	ms.randy->reset(time(NULL));
	ms.ReadHyper();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
