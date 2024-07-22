// Mutual_Inform_Meter.h: interface for the Mutual_Inform_Meter class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Histogram.h"
#include <stdio.h>
#include <math.h>
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class Mutual_Inform_Meter 
*/
/*! \ingroup Measurements 

\brief  Evaluates the mutual information of a sequence of Log Likelihood Ratio (LLR) vectors.

The \e istantaneous conditional entropy of a LLR vector \f$ \lambda_j(x) \f$ can be computed as
\f[
h_j ={\max_x}^*(\lambda_j(x))-\lambda_j(X_j)
\f]
where \f$ X_j \f$ is the transmitted symbol. The LLR vector has size \f$ M-1\f$, and by definition \f$ \lambda(0)=0 \f$.
LLR vectors are assumed to be in fixed point representation (integer values). The conversion to  real values
is done according the parameter #fact that is set in the method SetParameters()
\f[
\lambda = \lambda_I/fact
\f]

From the instantaneous conditional entropy,  the conditional entropy is computed as
\f[
H(X|Y)=\frac{1}{N}\sum_{i=0}^{N-1} h_j
\f]

The mutual information is computed as 
\f[
I(X;Y)=\log_2 M - H(X|Y)
\f]
All the above expressions are valid when the a-priori information of symbols is uniform.
*/
class Mutual_Inform_Meter  
{

public:
	
	Mutual_Inform_Meter();

	//! Set the main parameters.
	void SetParameters(
						const int M,			//!< Cardinality of input alphabet
						const double fact=8.	 //!< Factor for the quantization of LLR
						);

virtual ~Mutual_Inform_Meter();

	//! Compute istantaneous conditional entropy 
	/*! The istantaneous conditional entropy is computed as
	\f[
	h_i = {\max_x}^*(\lambda_i(x))-\lambda(X_i)\}
	\f]
	and all correspondent statistics are then updated. 
	If finite size analysis is not active (see method SetSizes()) then the istantaneous entropy
	is simply accumulated to compute a time average
	*/
	void Run(
		const int ntics,	//!< Number of symbols
		const int* X,		//!< index of transmitted symbols
		const int* lambda,	//!< Log-Likelihood ratios (one element is a vector of size M-1)
		const double* fact=0 //!< Optional sequence of factors for GEXIT computation (internal use)
		);

	//! Perform statistics for finite block sizes.
	/*! If this method is called, a statistic is performed also grouping blocks of \e N symbols.
	The histogram of the random variable
	\f[
	h_N(X|Y)=\frac{1}{N}\sum_{i=0}^{N-1} h_{jN+i}(X|Y)
	\f]
	is then evaluated. 
	This statistic can be used to compute upper and lower bounds to the frame error probability
	for finite block sizes.  
	*/
	void SetSizes(const int Nsizein,			//!< Number of considered sizes \e N
				const int* sizesin,				//!< Vector of sizes
				const int nrates,				//!< Number of rates considered 
				const double min,				//!< Minimum rate
				const double max				//!< Maximum rate
				);	
	
	//! Display results of statistic
	void Display(
		FILE* file=stdout,	//!< Output stream
		bool line=false		//!< flag to display results on a single line
		);

	//! Display the header
	void DisplayHeader(FILE* file=stdout);

	//! Reset all statistics
	void Reset();

	Histogram* hist;
	double GEXIT(){return gexit/((double)ic[0]*nsymbols);}
//	double Entropy(){return entropy/(fact*nsymbols);}
	double Entropy(){return entropy/((double)ic[0]*nsymbols);}
	double Dispersion(){
		return (V/((double)nsymbols)-(entropy/((double)nsymbols))*(entropy/((double)nsymbols)))/((double)ic[0]*ic[0]);}
	double MI(){return m-entropy/((double)ic[0]*nsymbols);}
	double gexit;
	double entropy;
	double V;

private:
	int* ic;
	int maxtab;
	int nsym;
	double m;	// Number of bits per symbol

	double fact;
	int minocc;
	__int64 nsymbols;
	double* totinf;
	int Nsize;	
	const int* sizes;
};
