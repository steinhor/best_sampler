#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "pratt_sampler/resonances.h"
Crandy *CresInfo::randy=NULL;
int CresInfo::NSPECTRAL=100;
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

void CresInfo::DecayGetResInfoPtr(int &nbodies,array<CresInfo *,5> &daughterresinfo){
	double r,bsum;
	int ibody,ibranch;
	CbranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	r=randy->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			printf("FATAL: In DecayGetResInfo: bsum too large, = %g\n",bsum);
			exit(1);
		}
	}while(bsum<r);
	nbodies=bptr->resinfo.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfo[ibody];
	}

}

bool CresInfo::CheckForDaughters(int pidcheck){
	//checks to see if any decay daughters match pid, for pid=0, checking for charged parts
	int ibody,nbodies,ibranch;
	bool exists=false;
	CresInfo *daughter;
	CbranchInfo *bptr;
	CbranchList::iterator bpos;
	if(pidcheck!=0){
		if(pid==pidcheck){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfo.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfo[ibody];
					if(daughter->pid==pidcheck){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(pidcheck)){
						exists=true;
						return exists;
					}
				}
				bpos++;
			}while(bpos!=branchlist.end());
		}
	}
	else{
		if(charge!=0){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfo.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfo[ibody];
					if(daughter->charge!=0){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(pidcheck)){
						exists=true;
						return exists;
					}
				}
				++bpos;
			}while(bpos!=branchlist.end());
		}
	}
	return exists;
}

CbranchInfo::CbranchInfo(){
}

void CresInfo::PrintBranchInfo(){
	Print();
	if(decay){
		printf(" ------  branches -------\n");
		for(int ib=0;ib<branchlist.size();ib++){
			printf("%2d, branching=%5.3f, ",ib,branchlist[ib]->branching);
			for(int ir=0;ir<branchlist[ib]->resinfo.size();ir++)
				printf("%6d ",branchlist[ib]->resinfo[ir]->pid);
			printf("Gamma_i=%g\n",width*branchlist[ib]->branching);
		}
	}
}

bool CresInfo::CheckForNeutral(){
	bool neutral=true;
	if(charge!=0 || strange!=0 || baryon!=0)
		neutral=false;
	return neutral;
}

void CresInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",pid,mass,minmass,name.c_str());
	printf("Gamma=%g, Degen=%g, Decay=%d\n",width,degen,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

void CresList::ReadResInfo(){
	//Cmerge *merge;
	int motherpid,pid,decay;
	double bsum,netm,bmax;
	int ires,ichannel,ibody,nbodies;
	int netq,netb,nets;
	string name, filename;
	CresInfo *resinfo=NULL,*aresinfo=NULL,*temp=NULL;
	CbranchInfo *bptr=NULL,*firstbptr=NULL,*abptr=NULL,*firstabptr=NULL;
	FILE *resinfofile;
	FILE *decayinfofile;
	char cname[200];
	bool antip;

	//dummy variables for decay info in resinfofile
	double gisospin;
	int dummy_int;
	int decays_Npart[25];
	double decays_branchratio[25];
	int decays_part[25][5];

	//dummy variables for resinfo while reading decays
	char d_name[200];
	double d_mass, d_width, d_gspin;
	int d_baryon, d_strange, d_charm, d_bottom;
	double d_gisospin;
	int d_charge,d_pid,d_L;

	filename=parmap->getS("RESONANCES_INFO_FILE",string("../local/resinfo/pdg-SMASH.dat"));
	printf("will read res info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	if (resinfofile==NULL) {
		fprintf(stderr,"Can't open resinfofile\n");
		printf("Error %d \n", errno);
		exit(1);
	}

	ires=0;
	int n=0;
	resinfo=new CresInfo();
	while(fscanf(resinfofile," %d",&d_pid)!=EOF && n<1000) {
		resinfo->pid=d_pid;
		n++;

		//main reading
		fscanf(resinfofile, " %s %lf %lf %lf %d %d %d %d %lf %d %d", cname,&d_mass,&d_width,&d_gspin,&d_baryon,&d_strange,&d_charm,&d_bottom,&gisospin,&d_charge,&decay);
		//using dummy variables so we can use the values again when creating the antiparticle
		resinfo->mass=d_mass;
		resinfo->width=d_width;
		resinfo->degen=d_gspin;
		resinfo->baryon=d_baryon;
		resinfo->strange=d_strange;
		resinfo->charm=d_charm;
		resinfo->bottom=d_bottom;
		resinfo->charge=d_charge;
		if(decay==1 && d_width==0.00000)
			{resinfo->decay=false;}
		else
			{resinfo->decay=true;}
		//resinfo->decay=bool((decay-1)); //subtract 1 because if stable, entry will be 1
		cname[int(strlen(cname))-1]='\0';
		resinfo->name=cname;

		//decay reading
		for (int j=0; j<decay; j++) { //reads into dummy variables: will read for decay info later
			fscanf(resinfofile, " %d %d %lf %d %d %d %d %d %d", &dummy_int,&decays_Npart[j],&decays_branchratio[j],&decays_part[j][0],&decays_part[j][1],&decays_part[j][2],&decays_part[j][3],&decays_part[j][4],&d_L);
		}

		if(resinfo->pid!=22){ //copied from old pid
			resinfo->ires=ires;
			ires+=1;
		}
		//else
		//	resinfo->ires=NResonances-1;
		resinfo->branchlist.clear();

		if(!resinfo->decay && (decay-1)==1){
			printf("decay turned off\n");
			resinfo->Print();
		}

		resmap.insert(CresInfoPair(resinfo->pid,resinfo));
		massmap.insert(CresMassPair(resinfo->mass,resinfo));


		//antiparticle creation
		if(resinfo->baryon!=0) {
			resinfo=new CresInfo();
			n+=1;

			//get antiparticle values from dummy variables created by original particle
			resinfo->pid=-d_pid;
			resinfo->mass=d_mass;
			resinfo->width=d_width;
			resinfo->degen=d_gspin;
			resinfo->baryon=-d_baryon;
			resinfo->strange=-d_strange;
			resinfo->charm=-d_charm;
			resinfo->bottom=-d_bottom;
			resinfo->charge=-d_charge;
			//resinfo->decay=bool((decay-1)); //subtract 1 because if stable, entry will be 1
			if(decay==1 && d_width==0.00000)
				{resinfo->decay=false;}
			else
				{resinfo->decay=true;}
			cname[int(strlen(cname))-1]='\0';
			string s(cname);
			resinfo->name="Anti-"+s;

			resinfo->ires=ires;
			ires+=1;

			resinfo->branchlist.clear();

			if(!resinfo->decay && (decay-1)==1){
				printf("decay turned off\n");
				resinfo->Print();
			}

			resmap.insert(CresInfoPair(resinfo->pid,resinfo));
			massmap.insert(CresMassPair(resinfo->mass,resinfo));
		}
		resinfo=new CresInfo();
	}
	printf("NResonances:%d\n",n);
	fclose(resinfofile);

	//Note: uses same file as resinfo, but needs to be read a second time because decay code needs to search resinfo map for products
	filename=parmap->getS("RESONANCES_DECAYS_FILE",string("../local/resinfo/pdg-SMASH.dat"));
	printf("will read decay info from %s\n", filename.c_str());
	decayinfofile=fopen(filename.c_str(), "r");
	if (resinfofile==NULL) {
		fprintf(stderr,"Can't open decayinfofile\n");
		printf("Error %d \n", errno);
		exit(1);
	}

	while(fscanf(decayinfofile," %d",&motherpid)!=EOF) {

		fscanf(decayinfofile," %s %lf %lf %lf %d %d %d %d %lf %d", d_name,&d_mass,&d_width,&d_gspin,&d_baryon,&d_strange,&d_charm,&d_bottom,&d_gisospin,&d_charge);
		fscanf(decayinfofile, " %d", &decay);
		//decay is the number of channels
		resinfo=GetResInfoPtr(motherpid);

		if (resinfo->baryon!=0){ //check if antiparticle exists
			antip=true;
			aresinfo=GetResInfoPtr(-motherpid); //resinfo for anti-particle
		}
		else
			antip=false;

		bsum=0.0;
		bmax=0.0;

		for (ichannel=0; ichannel<decay; ichannel++) {
			bptr=new CbranchInfo();
			bptr->resinfo.clear();
			resinfo->branchlist.push_back(bptr);

			if (antip) {
				abptr=new CbranchInfo(); //equivalent branch for the anti-particle
				abptr->resinfo.clear();
				aresinfo->branchlist.push_back(bptr);
			}

			fscanf(decayinfofile, " %d %d %lf", &dummy_int,&nbodies,&bptr->branching);

			if (antip)
				abptr->branching=bptr->branching;


			netq=-resinfo->charge;
			netb=-resinfo->baryon;
			nets=-resinfo->strange;
			netm=0.0;

			for(ibody=0; ibody<nbodies; ibody++) {
				fscanf(decayinfofile, " %d", &pid);
				temp=GetResInfoPtr(pid);
				bptr->resinfo.push_back(temp);

				if (antip){
					if(temp->baryon==0 && temp->charge==0 && temp->strange==0)
						abptr->resinfo.push_back(temp);
					else
						abptr->resinfo.push_back(GetResInfoPtr(-pid));
				}

				netq+=bptr->resinfo[ibody]->charge;
				netb+=bptr->resinfo[ibody]->baryon;
				nets+=bptr->resinfo[ibody]->strange;
				netm+=bptr->resinfo[ibody]->mass;
			}

			for (ibody=nbodies; ibody<5; ibody++) { //get rid of zeroes after real products
				fscanf(resinfofile, " %d", &dummy_int);
			}
			fscanf(resinfofile,"%d",&d_L); // angular momentum of decay
			bptr->L=d_L;

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
					if (antip) {
						//branching for antiparticle matches branching for particle
						firstabptr=aresinfo->branchlist[0];
						aresinfo->branchlist[0]=abptr;
						aresinfo->branchlist[ichannel]=firstabptr;
					}
				}
			}

		}  //out of channel loops
	}
	fclose(decayinfofile);
}

void CresInfo::CalcMinMass(){
	if(decay){
		minmass=1.0E20;
		double mbranch;
		int ibranch,nbodies,ibody;
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
