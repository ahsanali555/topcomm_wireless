/*! \file 
\brief Collects functions used to compute the symmetric capacity of a finitie size constellation 
 */
/*! \ingroup Design 
@{ */
#pragma once

// Returns the zeros of the Hermite functions of order n for quadrature formula
const double* HermiteZeroes(int n //!< Order of quadrature formula
	);

// Returns the weight of the Hermite functions of order n for quadrature formula
const double* HermiteWeights(int n//!< Order of quadrature formula
	);

// Returns the quadrature sum of fx
double Hermite(int n,	  //!< Order of the 
			   const double* fx //!< sum of values of the function at the symmetric nodes (f(x)+f(-x))
			   );

void GaussianIntegralnodes(
	const double fact,	//1/N0
	const int nn,		// Number of nodes
	double* n
	);

void GaussianIntegralnodes1D(const double fact, const int nn, double * n);

double GaussianIntegral1D(double * fx, const int nn);

double GaussianIntegral(const double* fx,const int nn);

// Information 
double Information(const double* n, 
		  const int m, 
		  const double *c, 
		  const double fact,
		  const bool symm=false);

// Return the Entropy of Y averaged on transmitted symbols
double Information_Lattice(
			const double* n,	 //!< Noise vector	
			const double *lv,	 //!< 2-nd vector defining the 2D-lattice (the first is assumed to be (1,0))
			const int m,		 //!< log2 of Number of points
			const double *c,	 //!< Constellation set
			const double fact,
			const bool symm
	);

// Pragmatic Information 
double Information_Prag(const double* n, 
		  const int m, 
		  const double *c, 
		  const double fact,
		  const bool symm=false);

//! Return the symmetric capacity and pragmatic capacity of a two-dimensional constellation.
/** The symmetric capacity is computed assuming equiprobable transmitted symbols.
The capacity is obtained using Hermite-Gauss quadrature rules for integration.
Capacity of a mod-lattice channel can also be computed
*/
double Capacity(
		  const int m,				//!< Log 2 of the cardinality of the constellation
		  const double *c,			//!< Pointer to Constellation points 
		  const double fact,		//!< inverse of N0
		  const bool prag=false,	//!< Compute pragmatic capacity
		  const int nn=7,			//!< Degree for Hermite interpolation
		  const bool symm=false,	//!< Apply symmetry condition
		  const double *lv=0		//!< Vector for defining 2D lattice (0 means no lattice)
		  );

//! Return the symmetric mutual indormation of a constellation for a AWGN channel affected by phase noise
/** The symmetric capacity is computed assuming equiprobable transmitted symbols.
The capacity is obtained using Hermite-Gauss quadrature rules for integration.
*/
double Inf_phasenoise(const int m,
		   const double*X,
		   const double* n,
		   const double Kn,
		   const double phi,
		   const double Kf,
		   const bool symm=false);

double Inf_phasenoise_prag(const int m,
		   const double*X,
		   const double* n,
		   const double Kn,
		   const double phi,
		   const double Kf,
		   const bool symm=false);
double MI_Phase(
				const int m,	//!< Modulation efficiency (number of bits)
				const double *c,	//!< constellation points (complex)
				const double factn,	//!< Concentration of additive noise (1/N0)
				const double factt,	//!< Concentration of phase jitter (\f$1/\sigma_\theta\f$)
				const bool prag=0,	//!< pragmatic MI?
				const int nn=3,		//!< Number of points for quadrature rule
				const bool symm=false		//!< Apply symmetry condition
				);	


double Information1D(const double* n, 
		  const int m, 
		  const double *c, 
		  const double fact);

// Pragmatic Information 
double Information_Prag1D(const double* n, 
		  const int m, 
		  const double *c, 
		  const double fact);

//! Return the symmetric capacity and pragmatic capacity of a one-dimensional constellation.
/** The symmetric capacity is computed assuming equiprobable transmitted symbols.
The capacity is obtained using Hermite-Gauss quadrature rules for integration.
*/
double Capacity1D(
		  const int m,			//!< Log 2 of the cardinality of the constellation
		  const double *c,		//!< Pointer to Constellation points 
		  const double fact,	//!< inverse of N0
		  const bool prag=false, //!< Compute pragmatic capacity
		  const int nn=3		//!< Degree for Hermite interpolation
		  );

double Delta_Information_Prag(const double* n, 
		  const int m, 
		  const double *c, 
		  const double fact);


double EntropyBAWGN(const double fact,const int nn=7);



void Entropy_vec(
	const int m,		//!< Log2 of number of points
	const double *c,	//!< constellation
	const double invn0,	//!< 1/N0
	double *H,			//!< vector of entropies
	const int nn=7		//!< Number of points in Hermite quadrature formula
	);	

void Entropy_vec2(
	const int m,		//!< Log2 of number of points
	const double *c,	//!< constellation
	const double invn0,	//!< 1/N0
	double *H,			//!< vector of entropies
	const int nn=7		//!< Number of points in Hermite quadrature formula
	);	

void Prag_Entropy_vec1D(const int m, const double * c, const int N, const int * mi, const double fact, double * H, const int nn);
void Prag_Information_vec1D(const double * n, const int m, const double * c, const int N, const int * mi, const double fact, double * acc, double * llr);
void Prag_Information_vec(const double * n, const int m, const double * c, const int N, const int * mi, const double fact, double * acc, double * llr);
void Prag_Entropy_vec(
		const int m,		//!< Log2 of number of points
		const double *c,	//!< constellation
		const int N,		//!< Number of stages
		const int* mi ,		//!< Bit per stage
		const double invn0, //!< 1/N0
		double *H,			//!< vector of entropies
		const int nn=7);

double gradient(
				const double* n,	// Noise 
				const int M,		// Number of points
				const double *c,	// points
				const double* modc,	// powers
				const double fact,	// 1/2/sigma^2
				const double tt,	// t for constraints
				double* g,			// gradient
				int average=1,		// Type on constraint on power
				int prg=0
				);

void QAWGN_Channel(
	const int m,		// Number of points
	const double *c,	// Constellation
	const double fact,
	double* f,      //! Transition matrix
	double* fn
	);


/*! @}*/