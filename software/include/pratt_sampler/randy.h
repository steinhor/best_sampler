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

class Crandy{
public:
	Crandy(int iseed);
	int seed;
	double ran();
	double ran_gauss();
	void ran_gauss2(double &g1,double &g2);
	double ran_exp();
	void reset(int iseed);
	void generate_boltzmann(double mass,double T,FourVector &p);
	void generate_boltzmann_alt(double mass,double T,FourVector &p);
	int poisson();
	void set_poisson_mean(double mu);  // For Poisson Dist
	
	int get_n(double delprob); // assuming each dp, such that integral dp =delprob,
	// gives probability dp of success, n = random number of successes
	// should be same as generating n for poissonian with mu=delprob
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


#endif