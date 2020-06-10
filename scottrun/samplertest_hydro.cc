#include "msu_sampler/master.h"
#include <cstring>
using namespace std;
using namespace msu_sampler;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("hydroparameters.txt");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	long long int deltacount=0,rhocount=0;
	CpartList pl=CpartList(&parmap);
	ms.partlist=&pl;
	ms.randy->reset(time(NULL));
	
	ms.ReadHyper();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS_TOT",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		deltacount+=ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
			+ms.partlist->CountResonances(2214)+ms.partlist->CountResonances(-2214)
				+ms.partlist->CountResonances(2224)+ms.partlist->CountResonances(-2224)
					+ms.partlist->CountResonances(1114)+ms.partlist->CountResonances(-1114);
		rhocount+=ms.partlist->CountResonances(113)+ms.partlist->CountResonances(213)+ms.partlist->CountResonances(-213);
		if((10*(ievent+1))%nevents==0)
			printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	printf("Nrho=%g\n",double(rhocount)/double(nevents));
	printf("Ndelta=%g\n",double(deltacount)/double(nevents));

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
