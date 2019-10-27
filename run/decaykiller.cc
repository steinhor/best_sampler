#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
#include <iostream>
#include <sstream>

using namespace std;

int main(){
	FILE *input=fopen("resonances_pdg.dat","r");
	FILE *output=fopen("resonances_nodecays.dat","w");
	int baryon,strangeness,charge,idecay,pid,Gparity;
	double mass,width,spin;
	char name[100];
	
	for(int iheader=0;iheader<5;iheader++){
		fgets(name,100,input);
		fprintf(output,"%s",name);
	}


	while(!feof(input)){
		fscanf(input,"%d %lf %d %d %d %lf %d %d %lf",&pid,&mass,&charge,&baryon,&strangeness,&spin,&Gparity,&idecay,&width);
		idecay=0;
		width=0.0;
		fgets(name,100,input);
		fprintf(output,"%d %g %d %d %d %g %d %d %g %s",pid,mass,charge,baryon,strangeness,spin,Gparity,idecay,width,name);
	}
	fclose(output);
	fclose(input);
	return 0;
}


