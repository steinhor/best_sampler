#include "pratt_sampler/master.h"
using namespace std;
using namespace pratt_sampler;

int main(int argc, char *argv[]){
	CparameterMap parmap;
	parmap.ReadParsFromFile("testparameters.txt");
	CmasterSampler ms(&parmap);
	int pid;
	printf("Enter pid: ");
	scanf("%d",&pid);
	CresInfo *resinfo=ms.reslist->GetResInfoPtr(pid);
	resinfo->PrintBranchInfo();
	resinfo->PrintSpectralFunction();
	return 0;
}
