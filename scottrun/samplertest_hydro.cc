#include "msu_sampler/master.h"
#include <cstring>
using namespace std;
using namespace msu_sampler;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("hydroparameters.txt");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	long long int count=0,anticount=0;
	CpartList pl=CpartList(&parmap);
	ms.partlist=&pl;
	ms.randy->reset(time(NULL));
	ms.ReadHyper();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS_TOT",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		count+=ms.partlist->CountResonances(2212)+ms.partlist->CountResonances(2112);
		anticount+=ms.partlist->CountResonances(-2212)+ms.partlist->CountResonances(-2112);
		if((ievent+1)%10==0)
			printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	printf("BBAR/B=%g\n",double(anticount)/double(count));

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
