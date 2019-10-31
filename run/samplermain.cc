#include "master.h"
#include "test.h"
#include <cstring>
using namespace std;

int main(int argc, char *argv[]){
	int nprint=100,iprint=0;

	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters.dat");

	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);

	ms.ReadHyper2D();

	//bool readout=false;
	//test::n_dummy(readout); //test function

	if (argc>1) {
		for (int i=1;i<argc;i++) {
			if (!strncmp(argv[i],"n_hyperlist",10)) test::n_hyperlist(ms,100);
			else if (!strncmp(argv[i],"n_dummy",10)) test::n_dummy(true);
			else if (!strncmp(argv[i],"bose_terms",10)) test::bose_terms();
			else if (!strncmp(argv[i],"test_viscous",10)) test::test_viscous();
			else if (!strncmp(argv[i],"newton_method_hyperlist",10)) test::newton_method_hyperlist(ms);
			else if (!strncmp(argv[i],"newton_method_dummy",10)) test::newton_method_dummy();
			else if (!strncmp(argv[i],"eps_rho_derivs",10)) test::eps_rho_derivs();
		}
	}

	/*****************************************************************************************
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	for(int ievent=0;ievent<nevents;ievent++){
		nparts+=ms.MakeEvent();
		printf("ievent=%d nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
	}
	*****************************************************************************************/

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
