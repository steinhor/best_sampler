#ifndef __HYPER_H__
#define __HYPER_H__
#include "classdefs.h"
#include "parametermap.h"
#include "randy.h"
#include "eos.h"
#include "misc.h"
using namespace std;

namespace msu_sampler {
	class Csampler;
	class CbulkQuantities;

	// ------------------------
	// Information about hyper-surface element
	//
	// ------------------------
	class Chyper{
	public:
		Chyper();
		//~Chyper();
		double T0; // T0 is read-in temperature, samplerptr->T is used to make parts;
		//Csampler *sampler;
		int firstcall;
		double sigma,rhoB,rhoI,rhoS;
		double muB,muI,muS,nhadrons,Rshear,Rbulk,RTbulk,epsilon,P,dedT;
		FourVector dOmega; // hyper volume
		double udotdOmega; // dot product of u and dOmega
		FourVector u;  // collective velocity
		FourVector r; // position
		FourVector qmu; // not used, but read in from Chun's hydro
		double tau,eta; // Bjorken tau and spatial rapidity
		FourTensor pitilde; // shear tensor
		double biggestpitilde; // largest eigenvalue
		void CalcBiggestpitilde();
		double PItilde; // bulk tensor correction
		double GetEntropyDensity();
		void SetSampler(Csampler *samplerptr);
		//double pixx,pixy,pixz,piyy,piyz,pizz;
		void FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz);
		void Print2D();  // prints out info for 2D element (Bjorken symm)
		void Print();    // prints out info for 3D element
		bool Rvisc_calculated;
		bool epsilon_calculated;
		Csampler *sampler;
	};
}

#endif
