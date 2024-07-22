#pragma once
/*! \file
\brief Collects functions used to optimize constellation sets he symmetric capacity of a finitie size constellation
*/
/*! \ingroup Design
@{ */
void Superposition_Modulation(const int n, const double *b, double *X, double& peak, double& ave);
void Superposition_Modulation(const int nd, const int n, const double *b, double *X, double& peak, double& ave);

void Quadrant_Modulation(const int n, const int nd, const double *X, double *XX, double& peak, double& ave);

void Generic_Superposition_Modulation(
	const int n,		//!< Number of services
	const int* bits,	//!< Bits in each service
	const double *b,	//!< Constellations in each service
	double *X,			//!< Output constellation
	double& peak,		//!< Peak power
	double& ave			//!< Average power
	);
void Optimize_Superposition_Modulation(
	const int m,			//! Number of modulation bits
	const bool optpeak,		//!< Flag to denote peak or average signal power optimization 
	const double invn0,		//!< SNIR (peak or average depending on previous flag) 
	double* &bbest,				//! Coefficients of best SP modulation [Output]
	double &MIbest,				//! The correspondent mutual information [Output]
	const int nmax=10000,		//! Number of SA steps  [optional]
	const double expt=1.,		//!< SA Cooling parameter		[optional]
	const double D0  =2.,		//!< SA initial deplacement	[optional]
	const double T0  =0.1		//!< SA initial temperature [optional]
	);




/* Program for the optimization of a finite size Modulation for a AWGN channel.
The program search for a constellation achieving the desired rate vector
*/
void Optimize_Modulation(
	const int m,			//!< Number of modulation bits
	const bool optpeak,		//!< Flag to denote peak or average signal power constraint 
	const double invn0,		//!< SNIR (peak or average depending on previous flag) 
	const double* XX,		//!< [in] Initial set of constellation points 
	double* Xbest,			//!< [out] Constellation points of best modulation 
	double &MIbest,			//!<[out] The correspondent mutual information [Output]
	const int nmax=10000,		//!<[in] Number of SA steps  
	const int nn=5,				//!<[in] 
	const bool quadrant=true,	//!< [in] apply quadrant symmetry
	const double expt=-0.4,		//!< [in] [SA Cooling parameter (default to polynomial cooling)
	const double D0  =2.,		//!< [in] SA initial deplacement	[optional]
	const double T0  =0.1,		//!< [in] SA initial temperature [optional]
	const double* lv=0   		//!< [in] Vector for lattice channel	
	);



/* Program for the optimization of a Hierarchical Modulation for a Broadcast Channel with service superposition
All services are superimposed  and a SIC receiver is assumed at the receiver.
Services are characterized by their (normalized) rates and maximal attenuations.
The program search for a constellation achieving the desired rate vector.
*/

int Optimize_HM(
	const int N,				//!< [in] Number of services
	const int mTOT,				//!< [in] Total number of bits in constellation
	const bool superp,			//!< [in] Superimposed (SP) constellation?
	const double Gamma,			//!< [in] Inverse of Noise power 
	const double GammaI,		//!< [in] Interference to Noise power	
	const double*  delta,		//!< [in] Difference of user attenuations in dB [a1/ax]
	const double*  rho,			//!< [in] Target Ratio of service rates [Rx/R1]
	const double* XX,			//!< [in] Initial Constellation (or SP coefficients) for optimization (XX =0 ->Random)	
	const int optpeak,			//!< [in] peak or average power constraint?
	const int optprag,			//!< [in] pragmatic constellation optimization?
	double &betabest,			//!< [out] best value of beta 
	double &betastd,			//!< [out] mean standard deviation w.r.t. target
	double* Rbest,				//!< [out] Rate of each service
	double* Xbest,				//!< [out] Best constellation (or SP coefficients)
	int* bits,					//!< [out] best bit allocation
	const int nmax=10000,		//!< [in] Number of trials Simulated Anneling		
	const double D0=1.,			//!< [in] Initial deplacement SA
	const double T0=0.1,		//!< [in] Initial temperature SA
	const double expt=-0.5,		//!< [in] type of cooling (expt<0 poynomial, 0<expt<1 exponential , expt=1 logarithmic)
	const bool quadrant=true,	//!< [in] Applies the quadrant condition on points
	const int nstarts=20,		//!< [in] Number of cold starts of SA
	const int lstart=200,		//!< [in] Length of cold starts SA		
	const int nn=5				//!< [in] Number of nodes for Gaussian quadrature formula
	);


/* Program for the optimization of a the Modulation for a two stage encoding
embedding a hard decoding stage.
The sum rate is maximized with the constraint of having the capacity of the hard decoded bits above 
a threshold specified by the user
*/

int Optimize_MLH(
	const int mTOT,				//!< [in] Total number of bits in constellation
	const int nbH,				//!< [in] Number of hard encoded bits
	const bool superp,			//!< [in] Superimposed (SP) constellation?
	const double Gamma,			//!< [in] Inverse of Noise power 
	const double* XX,			//!< [in] Initial Constellation (or SP coefficients) for optimization (XX =0 ->Random)	
	const int optpeak,			//!< [in] peak or average power constraint?
	const int optprag,			//!< [in] pragmatic constellation optimization?
	double &sumbest,			//!< [out] best value of sum rate 
	double* Rbest,				//!< [out] Rate of each service
	double* Xbest,				//!< [out] Best constellation (or SP coefficients)
	const int nmax=10000,		//!< [in] Number of trials Simulated Anneling		
	const double D0=1.,			//!< [in] Initial deplacement SA
	const double T0=0.1,		//!< [in] Initial temperature SA
	const double expt=-0.5,		//!< [in] type of cooling (expt<0 poynomial, 0<expt<1 exponential , expt=1 logarithmic)
	const bool quadrant=true,	//!< [in] Applies the quadrant condition on points
	const int nstarts=20,		//!< [in] Number of cold starts of SA
	const int lstart=200,		//!< [in] Length of cold starts SA		
	const int nn=5				//!< [in] Number of nodes for Gaussian quadrature formula
	);

/* Program for the optimization of the sum rate of a multistage encoder.
the Mutual information of each stage can be either the real one the pragmatic or the hard

*/

int Optimize_ML(
	const int N,				//!< [in] Number of services
	const int* bits,				//!< [in] Bits in each service
	const int* MItype,			//!< [in] Mutual information type to be optimized in each service
	const bool superp,			//!< [in] Superimposed (SP) constellation?
	const double Gamma,			//!< [in] Inverse of Noise power 
	const double* XX,			//!< [in] Initial Constellation (or SP coefficients) for optimization (XX =0 ->Random)	
	const int optpeak,			//!< [in] peak or average power constraint?
	double &sumbest,			//!< [out] best value of sum rate 
	double* Rbest,				//!< [out] Rate of each service
	double* Xbest,				//!< [out] Best constellation (or SP coefficients)
	const int nmax=10000,		//!< [in] Number of trials Simulated Anneling		
	const double D0=1.,			//!< [in] Initial deplacement SA
	const double T0=0.1f,		//!< [in] Initial temperature SA
	const double expt=-0.5,		//!< [in] type of cooling (expt<0 poynomial, 0<expt<1 exponential , expt=1 logarithmic)
	const bool quadrant=true,	//!< [in] Applies the quadrant condition on points
	const int nstarts=20,		//!< [in] Number of cold starts of SA
	const int lstart=200,		//!< [in] Length of cold starts SA		
	const int nn=5				//!< [in] Number of nodes for Gaussian quadrature formula
	);

int Optimize_MLSP(
	const int N,				//!< [in] Number of services
	const int* bits,				//!< [in] Bits in each service
	const int* MItype,			//!< [in] Mutual information type to be optimized in each service
	const double Gamma,			//!< [in] Inverse of Noise power 
	const double* XX,			//!< [in] Initial Constellation (or SP coefficients) for optimization (XX =0 ->Random)	
	const int optpeak,			//!< [in] peak or average power constraint?
	double &sumbest,			//!< [out] best value of sum rate 
	double* Rbest,				//!< [out] Rate of each service
	double* Xbest,				//!< [out] Best constellation (or SP coefficients)
	const int nmax = 100000,		//!< [in] Number of trials Simulated Anneling		
	const double D0 = 1.,			//!< [in] Initial deplacement SA
	const double T0 = 0.1f,		//!< [in] Initial temperature SA
	const double expt = -0.5,		//!< [in] type of cooling (expt<0 poynomial, 0<expt<1 exponential , expt=1 logarithmic)
	const int nstarts = 20,		//!< [in] Number of cold starts of SA
	const int lstart = 200,		//!< [in] Length of cold starts SA		
	const int nn = 5				//!< [in] Number of nodes for Gaussian quadrature formula
	);



double Optimize_Ortho(
	const double Nser,
	const double Gamma,
	const double GammaI,  
	const double*  delta, 
	const double*  rho,
	const int np,
	const double* C,
	double &betaortho
	);

int Optimize_ML1D(
	const int N, 
	const int * bits, 
	const int * MItype, 
	const bool superp, 
	const double invn0db, 
	const double * XX, 
	const int optpeak, 
	double & sumbest, 
	double * Rbest, 
	double * Xbest, 
	const int nmax, 
	const double D0, 
	const double T0, 
	const double expt, 
	const bool quadrant,  
	const int nstarts, 
	const int lstart, 
	const int nn);

/*! @}*/