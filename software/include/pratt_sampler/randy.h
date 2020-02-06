#ifndef __RANDY_H__
#define __RANDY_H__
#include <cstdlib>
#include <cmath>
#include <array>
#include <random>
#include "classdefs.h"

// ----------------------
// random number generation
// includes boltzmann generatione, e^{-E/T}, for parts of mass m
// thresholds are randomly placed points so that probability of 'netprob' passing any threshold in region epsilon->0 is epsilon ,
// regardless of past history.
// -----------------------

using namespace std;

namespace pratt_sampler {
	class Crandy{
	public:
		Crandy(int iseed);
		int seed;
		double ran();
		double ran_exp();
		void reset(int iseed);
		void generate_boltzmann(double mass,double T,FourVector &p);
		void generate_boltzmann_alt(double mass,double T,FourVector &p);
		bool test_threshold(double delprob); // if(netprob+delprob<threshold) increment netprob by delprob and return false,
		// else return true (used to see whether large sum of steps can be treated as one step)
		void increase_threshold();
		void increment_netprob(double delN);

		double netprob,threshold; // used for get_nsuccess
	private:
		std::mt19937 mt;
		std::uniform_real_distribution<double> ranu;
		std::normal_distribution<double> rang;
		std::poisson_distribution<int> ranp;
		//double netprob,threshold; // used for get_nsuccess
	};
}

#endif
