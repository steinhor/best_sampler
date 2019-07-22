#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "resonances.h"

Crandy *CresInfo::randy=NULL;

CresList::CresList(){
}

CresList::~CresList(){
	CresInfo *resinfo;
	CresInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		resmap.erase(rpos);
		delete resinfo;
		rpos=resmap.begin();
	}
	massmap.clear();
}

CresList::CresList(CparameterMap* parmap_in){
	parmap=parmap_in;
	//RESONANCE_DECAYS=parmap->getB("RESONANCE_DECAYS",true);
	ReadResInfo();
}

CresInfo::CresInfo(){
	minmass=0.0;
	branchlist.clear();
}

/*
Cspectralfunc::Cspectralfunc(){
}

Cspectralfunc_SMASH::Cspectralfunc_SMASH() : Cspectralfunc(){
}

double Cspectralfunc::rho(CresInfo *resinfo, double E, double kr)
{
    //Dummy
    return 0.0;
}

double Cspectralfunc_SMASH::rho(CresInfo *resinfo, double E, double kr){

    double rho,gamma,k;
    double resmass=resinfo->mass;
    double width=resinfo->width;
    double minmass=resinfo->minmass;
    double m1=resinfo->branchlist[0]->resinfo[0]->mass;
    double m2=resinfo->branchlist[0]->resinfo[1]->mass;
    double spin_deg=resinfo->spin;

    if(resinfo->branchlist[0]->resinfo[0]->decay==true || resinfo->branchlist[0]->resinfo[1]->decay==true)
    {
        double ma,mb,ma1,ma2,ma_pole,ma_0,ma_min,sum_ma,na,ma_gamma,ma_width;
        double form_lambda,ma_kr,ma_k,ma_rho,ma_rho0,suma,rho_width,rho_width_0,spectsum,spectsum0,ma_kra,ma_ka,s0;

        if(resinfo->branchlist[0]->resinfo[0]->decay==true)
        {   ma_min=resinfo->branchlist[0]->resinfo[0]->minmass;
            ma_pole=resinfo->branchlist[0]->resinfo[0]->mass;
            mb=resinfo->branchlist[0]->resinfo[1]->mass;
            ma_width=resinfo->branchlist[0]->resinfo[0]->width;
            ma1=resinfo->branchlist[0]->resinfo[0]->branchlist[0]->resinfo[0]->mass;
            ma2=resinfo->branchlist[0]->resinfo[0]->branchlist[0]->resinfo[1]->mass;
            if(m1==776 && m2==138) { form_lambda=0.8; }
            else if(resinfo->branchlist[0]->resinfo[1]->decay) { form_lambda=0.6; }
            else if(resinfo->branchlist[0]->resinfo[0]->baryon==0) { form_lambda=1.6; }
            else {form_lambda=2.0;}
        }
        if(resinfo->branchlist[0]->resinfo[1]->decay==true)
        {   ma_min=resinfo->branchlist[0]->resinfo[1]->minmass;
            ma_pole=resinfo->branchlist[0]->resinfo[1]->mass;
            mb=resinfo->branchlist[0]->resinfo[0]->mass;
            ma_width=resinfo->branchlist[0]->resinfo[1]->width;
            ma1=resinfo->branchlist[0]->resinfo[1]->branchlist[0]->resinfo[0]->mass;
            ma2=resinfo->branchlist[0]->resinfo[1]->branchlist[0]->resinfo[1]->mass;
            if(m1==776 && m2==138) { form_lambda=0.8; }
            else if(resinfo->branchlist[0]->resinfo[0]->decay) { form_lambda=0.6; }
            else if(resinfo->branchlist[0]->resinfo[1]->baryon==0) { form_lambda=1.6; }
            else {form_lambda=2.0;}
        }

        ma_kr=sqrt(abs(pow((ma_pole*ma_pole-ma1*ma1-ma2*ma2),2.0)-4.0*ma1*ma1*ma2*ma2))/(2.0*ma_pole);
        suma=0.0;
        int Na=100;

        for(int na=0;na<Na;na++)
        {
            double sum_ma=(na+0.5)/Na;
            ma_0 = 0.5*width*tan(PI*(sum_ma - .5));
            ma = ma_0+ma_pole;

            if(ma>=ma_min && ma<=E)
            {
                ma_k=sqrt(abs(pow((ma*ma-ma1*ma1-ma2*ma2),2.0)-(4.0*ma1*ma1*ma2*ma2)))/(2.0*ma);
                ma_gamma=ma_width*(ma_pole/ma)*((ma_k*ma_k*ma_k)/(ma_kr*ma_kr*ma_kr))*((ma_kr*ma_kr+HBARC*HBARC)/(ma_k*ma_k+HBARC*HBARC));
                ma_rho=(2.0)/(ma_width*PI)*0.25*ma_gamma*ma_gamma/((0.25*ma_gamma*ma_gamma)+(ma_pole-ma)*(ma_pole-ma));
                ma_rho0 = (1/PI)*(ma_width/2.0)/(0.25*ma_width*ma_width+ma_0*ma_0);
                ma_kra=sqrt(abs(pow((resmass*resmass-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*resmass);
                ma_ka=sqrt(abs(pow((E*E-ma*ma-mb*mb),2.0)-(4.0*ma*ma*mb*mb)))/(2.0*E);
                s0=ma+mb;
                rho_width=(ma_ka*ma_ka*ma_ka)/(E*(ma_ka*ma_ka+HBARC*HBARC))*((pow(form_lambda,4.0)+0.25*pow((s0-resmass*resmass),2.0))/(pow(form_lambda,4.0)+pow((E*E-0.5*(s0+resmass*resmass)),2.0)));
                rho_width_0=(ma_kra*ma_kra*ma_kra)/(resmass*(ma_kra*ma_kra+HBARC*HBARC));

                suma+=ma_rho/ma_rho0;
                spectsum+=rho_width*ma_rho/ma_rho0;
                spectsum0+=rho_width_0*ma_rho/ma_rho0;
            }
        }
        //if(ma>E) continue;
        double avg_weight_ma=suma/Na;
        double normal_ma=1.0/avg_weight_ma;
        double spect=normal_ma*spectsum/Na;
        double spect0=normal_ma*spectsum0/Na;
        gamma=width*spect/spect0;
    }

    else{
        k=sqrt(abs(pow((E*E-m1*m1-m2*m2),2.0)-(4.0*m1*m1*m2*m2)))/(2.0*E);
        if(spin_deg<1.001)
        {gamma=width*(resmass/E)*(k/kr); }
        else {gamma=width*(resmass/E)*((k*k*k)/(kr*kr*kr))*((kr*kr+HBARC*HBARC)/(k*k+HBARC*HBARC));}
    }

    rho=(2.0)/(width*PI)*0.25*gamma*gamma/((0.25*gamma*gamma)+(resmass-E)*(resmass-E));
    return rho;
}
*/

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

bool CresInfo::CheckForDaughters(int codecheck){
	//checks to see if any decay daughters match code, for code=0, checking for charged parts
	int ibody,nbodies,ibranch;
	bool exists=false;
	CresInfo *daughter;
	CbranchInfo *bptr;
	CbranchList::iterator bpos;
	if(codecheck!=0){
		if(code==codecheck){
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
					if(daughter->code==codecheck){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
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
					if(daughter->CheckForDaughters(codecheck)){
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

void CresInfo::DecayGetResInfoPtr_minmass(int &nbodies,array<CresInfo *,5> &daughterresinfo){
	nbodies=bptr_minmass->resinfo.size();
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr_minmass->resinfo[ibody];
	}
}

CbranchInfo::CbranchInfo(){
}

bool CresInfo::CheckForNeutral(){
	bool neutral=true;
	if(charge!=0 || strange!=0 || baryon!=0)
		neutral=false;
	return neutral;
}

void CresInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",code,mass,minmass,name.c_str());
	printf("Gamma=%g, Spin=%g, Decay=%d\n",width,spin,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

void CresList::ReadResInfo(){
	//Cmerge *merge;
	int mothercode,code,decay,strange,charge,baryon,NResonances;
	double mass,mothermass,spin,width,bsum,netm,qR2,bmax;
	int iires,ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,length, LDecay, i_inel;
	int netq,netb,nets, netg, G_Parity;
	string name, filename;
	CresInfo *resinfo=NULL,*oldresinfo=NULL, *resinfo_1 = NULL, *aresinfo=NULL, *temp=NULL;
	CbranchInfo *bptr=NULL,*oldbptr=NULL,*firstbptr=NULL,*abptr=NULL, *firstabptr=NULL;
	FILE *resinfofile;
	FILE *decayinfofile;
	char dummy[200],cname[200];
	bool antip;

	//dummy variables for decay info in resinfofile
	double gisospin;
	double gspin;
	int dummy_int;
	int decays_Npart[20];
	double decays_branchratio[20];
	int decays_part[20][5];

	//dummy variables for resinfo while reading decays
	char d_name[200];
	double d_mass, d_width, d_gspin;
	int d_baryon, d_strange, d_charm, d_bottom;
	double d_gisospin;
	int d_charge,d_code;

	filename=parmap->getS("RESONANCES_INFO_FILE",string("../local/resinfo/pdg-SMASH.dat"));
	printf("will read res info from %s\n",filename.c_str());

	/************************************************************************************************
	#OLD RESINFO READING CODE:
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	printf("NResonances=%d\n",NResonances);
	ires=0;
	for(iires=0;iires<NResonances;iires++){
		resinfo=new CresInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfo->code,&resinfo->mass,&resinfo->charge,&resinfo->baryon, &resinfo->strange,&resinfo->spin,&resinfo->G_Parity,&decay,&resinfo->width);
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfo->name=cname;
		resinfo->decay=bool(decay);
		//if(!RESONANCE_DECAYS)
		//	resinfo->decay=false;
		if(resinfo->code!=22){
			resinfo->ires=ires;
			ires+=1;
		}
		else
			resinfo->ires=NResonances-1;
		resinfo->branchlist.clear();
		if(!resinfo->decay && decay==1){
			printf("decay turned off\n");
			resinfo->Print();
		}
		resmap.insert(CresInfoPair(resinfo->code,resinfo));
		massmap.insert(CresMassPair(resinfo->mass,resinfo));
	}
	fclose(resinfofile);
	************************************************************************************************/

	resinfofile=fopen(filename.c_str(),"r");
	if (resinfofile==NULL) {
		fprintf(stderr,"Can't open resinfofile\n");
		printf("Error %d \n", errno);
		exit(1);
	}

	ires=0;
	int n=0;
	resinfo=new CresInfo();
	while(fscanf(resinfofile," %d",&d_code)!=EOF && n<1000) {
		resinfo->code=d_code;
		n++;

		//main reading
		fscanf(resinfofile, " %s %lf %lf %lf %d %d %d %d %lf %d %d", cname,&d_mass,&d_width,&d_gspin,&d_baryon,&d_strange,&d_charm,&d_bottom,&gisospin,&d_charge,&decay);
		//using dummy variables so we can use the values again when creating the antiparticle
		resinfo->mass=d_mass;
		resinfo->width=d_width;
		resinfo->spin=d_gspin;
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
			fscanf(resinfofile, " %d %d %lf %d %d %d %d %d", &dummy_int,&decays_Npart[j],&decays_branchratio[j],&decays_part[j][0],&decays_part[j][1],&decays_part[j][2],&decays_part[j][3],&decays_part[j][4]);
		}

		if(resinfo->code!=22){ //copied from old code
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

		resmap.insert(CresInfoPair(resinfo->code,resinfo));
		massmap.insert(CresMassPair(resinfo->mass,resinfo));


		//antiparticle creation
		if(resinfo->baryon!=0) {
			resinfo=new CresInfo();
			n++;

			//get antiparticle values from dummy variables created by original particle
			resinfo->code=-d_code;
			resinfo->mass=d_mass;
			resinfo->width=d_width;
			resinfo->spin=d_gspin;
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

			resmap.insert(CresInfoPair(resinfo->code,resinfo));
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

	while(fscanf(decayinfofile," %d",&mothercode)!=EOF) {

		fscanf(decayinfofile," %s %lf %lf %lf %d %d %d %d %lf %d", d_name,&d_mass,&d_width,&d_gspin,&d_baryon,&d_strange,&d_charm,&d_bottom,&d_gisospin,&d_charge);
		fscanf(decayinfofile, " %d", &decay);
		//decay is the number of channels

		resinfo=GetResInfoPtr(mothercode);
		resinfo->minmass=1.0E10;

		if (resinfo->baryon!=0) { //check if antiparticle exists
			antip=true;
			aresinfo=GetResInfoPtr(-mothercode); //resinfo for anti-particle
			aresinfo->minmass=1.0E10;
		} else {antip=false;}

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

			if (antip) {
				abptr->branching=bptr->branching;
			}

			netq=-resinfo->charge;
			netb=-resinfo->baryon;
			nets=-resinfo->strange;
			netm=0.0;

			for(ibody=0; ibody<nbodies; ibody++) {
				fscanf(decayinfofile, " %d", &code);
				temp=GetResInfoPtr(code);
				bptr->resinfo.push_back(temp);

				if (antip) {
					if(temp->baryon==0 && temp->charge==0 && temp->strange==0){
						abptr->resinfo.push_back(temp);
					}
					else {abptr->resinfo.push_back(GetResInfoPtr(-code));}
				}

				netq+=bptr->resinfo[ibody]->charge;
				netb+=bptr->resinfo[ibody]->baryon;
				nets+=bptr->resinfo[ibody]->strange;
				netm+=bptr->resinfo[ibody]->mass;
			}

			for (ibody=nbodies; ibody<5; ibody++) { //get rid of zeroes after real products
				fscanf(resinfofile, " %d", &dummy_int);
			}

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

			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfo->minmass){
				resinfo->minmass=netm;
				resinfo->bptr_minmass=bptr;

				if (antip) {
					aresinfo->minmass=netm;
					aresinfo->bptr_minmass=abptr;
				}
			}
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

		}

	}
	fclose(decayinfofile);
	/*********************************************************************************************************
	//OLD DECAY READING CODE:
	filename=parmap->getS("RESONANCES_DECAYS_FILE",string("resinfo/decays_pdg_weak.dat"));
	printf("will read decay info from %s\n",filename.c_str());
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfo=GetResInfoPtr(mothercode);
		resinfo->minmass=1.0E10;
		bsum=0.0;
		bmax=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CbranchInfo();
			bptr->resinfo.clear();
			resinfo->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfo->charge;
			netb=-resinfo->baryon;
			nets=-resinfo->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfo.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfo[ibody]->charge;
				netb+=bptr->resinfo[ibody]->baryon;
				nets+=bptr->resinfo[ibody]->strange;
				netm+=bptr->resinfo[ibody]->mass;
			}

			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				printf("Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				printf("MOTHER (ichannel=%d, nbodies=%d):\n",ichannel,nbodies);
				resinfo->Print();
				printf("DAUGHTERS:\n");
				for(ibody=0;ibody<nbodies;ibody++)
					bptr->resinfo[ibody]->Print();
				if(netq!=0 || netb!=0)
					exit(1);
			}
			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			//store two body decays only
			//Note: there was an open multi-line comment symbol here!!!!!!!!!
			if(nbodies==2){
				ires1=bptr->resinfo[0]->ires;
				ires2=bptr->resinfo[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new Cmerge(resinfo,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new Cmerge(resinfo,bptr->branching, LDecay);
				}

				if(resinfo->mass>bptr->resinfo[0]->mass+bptr->resinfo[1]->mass){
					qR2=Misc::triangle(resinfo->mass,
					bptr->resinfo[0]->mass,bptr->resinfo[1]->mass);
					SigmaMaxArray[ires1][ires2]+=
						(bptr->branching*4.0*PI*HBARC*HBARC)*(2.0*resinfo->spin+1.0)/
							((2.0*bptr->resinfo[0]->spin+1.0)*(2.0*bptr->resinfo[1]->spin+1.0)*qR2);
					if(ires2!=ires1)
						SigmaMaxArray[ires2][ires1]=SigmaMaxArray[ires1][ires2];

				}
			}
			//Note: there was a close multi-line comment symbol here!!!!!!!!!
			bsum+=bptr->branching;
			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfo->minmass){
				resinfo->minmass=netm;
				resinfo->bptr_minmass=bptr;
			}
			// switch places to make sure first branch has largest
			if(bptr->branching>bmax){
				bmax=bptr->branching>bmax;
				if(ichannel>0){
					firstbptr=resinfo->branchlist[0];
					resinfo->branchlist[0]=bptr;
					resinfo->branchlist[ichannel]=firstbptr;
				}
			}
		}
	}
	fclose(decayinfofile);
	*********************************************************************************************************/
}

CresInfo* CresList::GetResInfoPtr(int code){
	CresInfoMap::iterator rpos;
	rpos=resmap.find(code);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		printf("Warning GetResInfoPtr() can't find match for PID=%d\n",code);
		exit(1);
		return NULL;
	}
}

#endif
