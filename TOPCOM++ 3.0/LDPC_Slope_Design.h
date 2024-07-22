#include "LDPC_Encoder.h"
#pragma once
/*! \file
\brief Collects functions used to optimize degree distributions of LDPC codes with constrained complexity and slope
*/
/*! \ingroup Design
@{ */
double eval(const double * L, const double x);
double Eval(const double * L, const double x);
double Norm(double * L);
double Rk(const int k, const double * R);
void Regular_degree(const double rate, const double eta, double *& L, double *& R, const int maxdeg);
double Slope(const double eot, const double eta, const double rate, const int maxdeg);
double densityatslope(const double eot, const double slopet, const double rate, const int maxdeg);
bool accdegree(double *& Lacc, double * L, const double f=1.,const int xmin=1);
double Threshold(const double * L, const double * R, const double slope, double * eio, double * eoo, double * xo);
double ThresholdP(const double * L, const double * R, const double slope, const int P, double * eio, double * eoo, double * xo);
void GetCode1(const double slope, const double rate, int nterms, double *& R, double *& L);
double GetCodewithDensity(
	double *& L,			//!< Left degree profile (node perspective)
	double *& R,			//!< Right degree profile (node perspective)
	int & degstepbest,		//!< Delta of best solution
	const double rate,		//!< Target rate
	const double S,			//!< Target slope
	const double etat,		//!< Target density
	const int degstep=0,	//!< Chosen step of degree (0 for optimal)
	const int maxdeg=10000,	//!< Maximum variable degree
	const double eps = 0.01,	  //!< Step for increasing coefficient (negative for not normalizing to d)
	const double strongstab = 1., //!< Factor for stability condition
	const int P = 1000000        //!< Minimum degree of punctured bit
	);

//! Return the degree distribution of code with minimum density (but below some maximal value) and threshold above given value
double GetCodewithThreshold(
	double *& L,			//!< Left degree profile (node perspective)
	double *& R,			//!< Right degree profile (node perspective)
	int & degstepbest,		//!< Delta of best solution
	const double rate,		//!< Target rate
	const double S,			//!< Target slope
	const double alphat,		//!< Target code efficiency
	const double etamax,	//!< Maximum density
	const int degstep=0,	//!< Chosen step of degree (0 for optimal)
	const int maxdeg=10000,	//!< Maximum variable degree
	const double eps=0.01,	//!< Step for increasing coefficient (negative for not normalizing to d)
	const double strongstab = 0. //!< Factor for stability condition
	);


void DisplayEXIT(double * L, double * R, FILE * file);
void DisplayLREXIT(double* L, double* R, double ei, FILE* file);
void Display(double * L, double * R, FILE * file=stdout, bool online=false);

//! Returns a Good Linearly encodable LDPC with target K, N, lifting factor Z, Slope, and complexity
/*!
- Computes optimal left degree distribution with target complexity an slope
- Maximize girth fow a given K,N, and lifting factor Z
- Makes the code linearly encodable applying the richardson algorithm
*/

LDPC_Encoder* GetGoodLDPC(
	const int K,					
	const int N,
	const int Z,
	const double etat = 3.5,		//!< Target density 
	const double stabfact = 0.,	    //!< Stabilty factor (low error floors)
	const bool withP = false,		//!< Allows puncturing of some variable nodes
	const double S = 0.0,			//!< Target slope
	const bool PEG = true,			//!< Applies PEG
	double* th=0					//!< Optional pointer  to store code efficiency
	);


double EvalMN(
	const double * NU, // Multimomial (number of terms,number of variables, coeff and degrees)
	const double *x = 0,  // default to all one
	const int ix = 0,	  // derivative w.r.t. one variable
	const bool direct = true	  // (use x or 1-x)
	);

void MEDisplay(const double* L, FILE* ff = stdout);
double* ME_Fixed(const double * NU, const double * MU, const double x);
double ME_Threshold(const double * NU, const double * MU,double*x=0);
double GetMECodewithDensity(
	double *& NU,		//!< [out] Left degree of optimized Multi-edge distribution
	double *& MU,       //!< [out] Right degree of optimized Multi-edge distribution
	const double R,		//!< [in]  Rate of the designed code
	double &Rx,         //!< [in]  Rate of core code (Rx>R)
	const double Cmax=100,
	const double eps=0.01,
	const int maxdeg=100	//!< Maximum variable node degree 
	);
/*! @}*/
