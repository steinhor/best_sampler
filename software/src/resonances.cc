#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "pratt_sampler/resonances.h"
using namespace pratt_sampler;

Crandy *CresInfo::randy=NULL;
unsigned int CresInfo::NSPECTRAL=100;
string CresInfo::SFDIRNAME="../local/resinfo/spectralfunctions";

CresList::CresList(){
}

CresList::~CresList(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		resmap.erase(rpos);
		rpos=resmap.begin();
		delete resinfo;
	}
	resmap.clear();
	massmap.clear();
}

CresList::CresList(CparameterMap* parmap_in){
	parmap=parmap_in;
	CresInfo::NSPECTRAL=parmap->getI("SAMPLER_NSPECTRAL",100);
	CresInfo::SFDIRNAME=parmap->getS("SAMPLER_SFDIRNAME","../local/resinfo/spectralfunctions");
	//RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	ReadResInfo();
	CalcMinMasses();
	//CalcSpectralFunctions();
	ReadSpectralFunctions();
}

CresInfo::CresInfo(){
	minmass=1.0E20;
	SFcalculated=false;
	branchlist.clear();
}

void CresList::CalcMinMasses(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos;
	for(rpos=resmap.begin();rpos!=resmap.end();++rpos){
		resinfo=rpos->second;
		resinfo->CalcMinMass();
	}
}

CbranchInfo::CbranchInfo(){
}

void CresInfo::PrintBranchInfo(){
	Print();
	if(decay){
		printf(" ------  branches -------\n");
		for(unsigned int ib=0;ib<branchlist.size();ib++){
			printf("%2du, branching=%5.3f, ",ib,branchlist[ib]->branching);
			for(unsigned int ir=0;ir<branchlist[ib]->resinfo.size();ir++)
				printf("%6du ",branchlist[ib]->resinfo[ir]->pid);
			printf("Gamma_i=%g\n",width*branchlist[ib]->branching);
		}
	}
}

void CresInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",pid,mass,minmass,name.c_str());
	printf("Gamma=%g, Degen=%g, Decay=%d\n",width,degen,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

void CresList::ReadResInfo(){
	//Cmerge *merge;
	int motherpid,pid;
	double bsum,netm,bmax;
	unsigned int ires,ichannel,ibody,nbodies,n;
	int netq,netb,nets;
	string name, filename;
	CresInfo *resinfo=NULL,*aresinfo=NULL,*temp=NULL;
	CdecayInfo *decayinfo, *adecayinfo;
	CbranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	char cname[200];
	CdecayInfoMap decaymap;
	double gisospin;
	int dummy_int;
	CresInfoMap::iterator iter;

	filename=parmap->getS("RESONANCES_INFO_FILE",string("../software/resinfo/pdg-SMASH.dat"));
	printf("will read resonance info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	if (resinfofile==NULL) {
		fprintf(stderr,"Can't open resinfofile\n");
		printf("Error %d \n", errno);
		exit(1);
	}

	ires=0;
	n=0;
	while(fscanf(resinfofile," %d",&pid)!=EOF && n<1000) {
		resinfo=new CresInfo();
		decayinfo=new CdecayInfo();
		n++;
		resinfo->pid=pid;

		//main reading
		fscanf(resinfofile, " %s %lf %lf %lf %d %d %d %d %lf %d %d", cname,&resinfo->mass,&resinfo->width,&resinfo->degen,&resinfo->baryon,&resinfo->strange,&resinfo->charm,&resinfo->bottom,&gisospin,&resinfo->charge,&resinfo->nchannels);
		if(resinfo->nchannels==1 && resinfo->width==0.00000)
			{resinfo->decay=false;}
		else
			{resinfo->decay=true;}
		//resinfo->decay=bool((decay-1)); //subtract 1 because if stable, entry will be 1
		cname[int(strlen(cname))-1]='\0';
		resinfo->name=cname;

		//decay reading
		for (unsigned int j=0; j<resinfo->nchannels; j++) { //reads into map values: will access for decays when done creating resonances
			fscanf(resinfofile, " %d %d %lf %d %d %d %d %d %d", &dummy_int,&decayinfo->Nparts[j],&decayinfo->branchratio[j],&decayinfo->products[j][0],&decayinfo->products[j][1],&decayinfo->products[j][2],&decayinfo->products[j][3],&decayinfo->products[j][4],&decayinfo->d_L[j]);
		}

		if(resinfo->pid!=22){ //copied from old pid
			resinfo->ires=ires;
			ires+=1;
		}

		resinfo->branchlist.clear();

		if(!resinfo->decay && (resinfo->nchannels-1)==1){
			printf("decay turned off\n");
			resinfo->Print();
		}

		resmap.insert(CresInfoPair(resinfo->pid,resinfo));
		massmap.insert(CresMassPair(resinfo->mass,resinfo));
		decaymap.insert(CdecayInfoPair(resinfo->pid,decayinfo));

		//antiparticle creation
		if(resinfo->baryon!=0) {
			aresinfo=new CresInfo();
			adecayinfo=new CdecayInfo();
			n+=1;

			aresinfo->pid=-resinfo->pid;
			aresinfo->mass=resinfo->mass;
			aresinfo->width=resinfo->width;
			aresinfo->degen=resinfo->degen;
			aresinfo->baryon=-resinfo->baryon;
			aresinfo->strange=-resinfo->strange;
			aresinfo->charm=-resinfo->charm;
			aresinfo->bottom=-resinfo->bottom;
			aresinfo->charge=-resinfo->charge;
			aresinfo->decay=resinfo->decay;
			aresinfo->nchannels=resinfo->nchannels;
			cname[int(strlen(cname))-1]='\0';
			string s(cname);
			aresinfo->name="Anti-"+s;
			aresinfo->ires=ires;
			ires+=1;

			aresinfo->branchlist.clear();

			for (unsigned int j=0; j<resinfo->nchannels; j++) { //reads into map values: will access for decays when done creating resonances
				for (unsigned int i=0; i<5; i++) {
					pid=decayinfo->products[j][i];
					if (pid!=0) {
						temp=GetResInfoPtr(pid);
						if(temp->baryon==0 && temp->charge==0 && temp->strange==0){
							adecayinfo->products[j][i]=decayinfo->products[j][i];
						}
						else {
							adecayinfo->products[j][i]=-decayinfo->products[j][i];
						}
					}
					else adecayinfo->products[j][i]=0;
				}
				adecayinfo->Nparts[j]=decayinfo->Nparts[j];
				adecayinfo->branchratio[j]=decayinfo->branchratio[j];
				adecayinfo->d_L[j]=decayinfo->d_L[j];
			}
			resmap.insert(CresInfoPair(aresinfo->pid,aresinfo));
			massmap.insert(CresMassPair(aresinfo->mass,aresinfo));
			decaymap.insert(CdecayInfoPair(aresinfo->pid,adecayinfo));
		}
	}
	printf("NResonances:%d\n",n);
	fclose(resinfofile);

	//now, use the stored decay information to create branchlists
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		motherpid=iter->first;
		decayinfo=decaymap[motherpid];

		bsum=0.0;
		bmax=0.0;

		for (ichannel=0; ichannel<resinfo->nchannels; ichannel++) {
			nbodies=decayinfo->Nparts[ichannel];
			bptr=new CbranchInfo();
			bptr->resinfo.clear();
			resinfo->branchlist.push_back(bptr);
			bptr->branching=decayinfo->branchratio[ichannel];

			netq=-resinfo->charge;
			netb=-resinfo->baryon;
			nets=-resinfo->strange;
			netm=0.0;

			for(ibody=0; ibody<nbodies; ibody++) {
				pid=decayinfo->products[ichannel][ibody];
				bptr->resinfo.push_back(GetResInfoPtr(pid));
				netq+=bptr->resinfo[ibody]->charge;
				netb+=bptr->resinfo[ibody]->baryon;
				nets+=bptr->resinfo[ibody]->strange;
				netm+=bptr->resinfo[ibody]->mass;
			}

			bptr->L=decayinfo->d_L[ichannel];

			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				printf("Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				resinfo->Print();
				printf("DAUGHTERS:\n");
				for(ibody=0;ibody<nbodies;ibody++)
					bptr->resinfo[ibody]->Print();
				if(netq!=0 || netb!=0)
					exit(1);
			}

			bsum+=bptr->branching;

			// switch places to make sure first branch has largest
			if(bptr->branching>bmax){
				bmax=bptr->branching;
				if(ichannel>0){
					firstbptr=resinfo->branchlist[0];
					resinfo->branchlist[0]=bptr;
					resinfo->branchlist[ichannel]=firstbptr;
				}
			}
		}  //out of channel loops
	}
	for(iter=resmap.begin();iter!=resmap.end();++iter){
		resinfo=iter->second;
		motherpid=iter->first;
		decayinfo=decaymap[motherpid];
		delete decayinfo;
	}
	decaymap.clear();
}

void CresInfo::CalcMinMass(){
	if(decay){
		minmass=1.0E20;
		double mbranch;
		unsigned int ibranch,nbodies,ibody;
		CbranchInfo *bptr;
		for(ibranch=0;ibranch<branchlist.size();ibranch++){
			bptr=branchlist[ibranch];
			mbranch=0.0;
			nbodies=bptr->resinfo.size();
			for(ibody=0;ibody<nbodies;ibody++){
				if(bptr->resinfo[ibody]->minmass>1.0E8){
					bptr->resinfo[ibody]->CalcMinMass();
				}
				mbranch+=bptr->resinfo[ibody]->minmass;
			}
			if(mbranch<minmass)
				minmass=mbranch;
		}
	}
	else
		minmass=mass;
}

CresInfo* CresList::GetResInfoPtr(int pid){
	CresInfoMap::iterator rpos;
	rpos=resmap.find(pid);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		printf("Warning GetResInfoPtr() can't find match for PID=%d\n",pid);
		exit(1);
		return NULL;
	}
}

#endif
