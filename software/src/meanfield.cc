#include "meanfield.h"

CmeanField::CmeanField(){
	//
}

CmeanField_Simple::CmeanField_Simple(CparameterMap *parmap) : CmeanField(){
}

double CmeanField_Simple::GetMass(CresInfo *resinfo,double sigma){
	if(resinfo->baryon!=0)
		return (sigma/93.0)*resinfo->mass;
	else
		return resinfo->mass;
}

double CmeanField::GetMass(CresInfo *resinfo,double sigma){
	// Dummy function
	return 0.0;
}
