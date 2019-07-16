#ifndef __HYPER_H__
#define __HYPER_H__
#include "classdefs.h"
#include "constants.h"
#include "parametermap.h"
#include "randy.h"
#include "eos.h"
#include "misc.h"

using namespace std;
class Csampler;
class CbulkQuantities;

class Chyper{
public:
	Chyper();
	int ihyp;
	double T,sigma,rhoB,rhoI,rhoS;
	double muB,muI,muS,nhadrons,lambda,epsilon;
	FourVector dOmega; // hyper volume
	double udotdOmega;
	FourVector u;  // collective velocity
	FourVector r; // position
	double tau,eta; // Bjorken tau and spatial rapidity
	double pitilde[4][4]; // shear tensor
	//double pixx,pixy,pixz,piyy,piyz,pizz;
	void FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz);
	void Print2D();  // prints out info for 2D element (Bjorken symm)
	void Print();    // prints out info for 3D element
};

#endif
