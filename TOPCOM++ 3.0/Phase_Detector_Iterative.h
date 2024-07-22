// Phase_Detector_Iterative.h: interface for the Phase_Detector_Iterative class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Modulator.h"
#include "Pilot.h"
#include "ctopcom.h"
#include "Include.h"
/*! \file 
\brief Declaration of class Phase_Detector_Iterative
*/
/*! \ingroup Synchronization
\brief Iterative Phase detector based on Colavolpe and Caire algorithm

The class  implements a block for the compensation of phase jitter based on an iterative
algorithm that exploits pilot symbols and a-priori information on the
transmitted symbols. The block can be used iteratively in connection
with an outer SISO channel decoder, and provides very good results
with low signal to noise ratio.

The algorithm is the ``recursive algorithm", based on a forward and backward recursion.

The received complex samples are assumed to be
\f[
r_k=m(c_k) e^{-j\theta_k}+\sigma n_k,
\f]
where \f$m\f$ is the modulator mapping that associates complex
constellation points to an index in \f$c_k\in\{0,\ldots,M-1\}\f$. The
underlying assumed phase process \f$\theta_k\f$ is a Wiener process
\f[
\theta_{k+1}= \theta_k +\sigma_\Delta t_k
\f]
and  \f$t_k\f$ and \f$n_k\f$ are white discrete Gaussian processes with unit
variance and zero mean. The model  is fully described by the
variances \f$\sigma_\Delta^2\f$ of the phase step and \f$\sigma\f$.

The block accept and provide soft-information about the data symbols
in the form of quantized LLR, vectors of \f$M-1\f$ dimensions defined as
\f[
\lambda_k(i) =\log\frac{p(c_k=i+1)}{p(c_k=0)} \;i=0,\ldots,M-2
\f]

The block is configured SetParameters() with a pointer to the
reference Modulator for the data, a pointer to a Pilot block
that specifies the positions and type of pilots, the variance of the
phase steps of Wiener process \f$\sigma_\Delta\f$,  the variance of the
additive white noise \f$\sigma_\Delta\f$, the frame length \f$K\f$, including
the pilots, and finally the factor used for quantization of LLR,
which is assumed by default to be 8.

A tic of the Run method process a whole block of data and pilots.
The number of data in a frame \f$K_d\f$ is internally computed by
extracting the pilots symbols from the whole frame \f$K\f$.

Input signals are the \f$K\f$  received samples affected by phase jitter
and the set of \f$K_d\f$ input LLR of the data symbols. The
a-priori knowledge of pilot symbols is exploited internally and the
user does not need to provide a-priori about the pilots. If no
a-priori information about the data symbols is available, the user
can specify the null pointer.

Output signals are the updated \f$K_d\f$ LLR of the transmitted symbols.

For an example of its use see e.g. the test program "test_iterative_phase.cpp".

  \author Guido Montorsi
*/

class Phase_Detector_Iterative  
{
public:
	Phase_Detector_Iterative();

	//! Destructor
	virtual ~Phase_Detector_Iterative();

	//! Set the parameters of the algorithm
	void SetParameters(
		const Modulator *mod,   //!< Reference to the used Modulator
		Pilot *pil,				//!< Reference to the pilot sequence
		double sigmaD,			//!< Variance of phase steps of Wiener process
		double sigma,			//!< Variance of additive noise
		const int K,			//!< Block length (number of symbols)
		const double fact=8.	//!< Factor for quantization of LLR
		);

	//! Simulation method
	void Run(const int tics,   //!< Number of processed samples
		const double *r,	   //!< Received samples, affected by phase jitter
		const int *llri,	   //!< LLR of constellation points
		int* llro			    //!< Updated LLR
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,		//!< Number of processed samples
		const cmplx* r,				//!< Complex Received samples, affected by phase jitter
		const int *llri,			//!< LLR of constellation points
		int* llro					//!< Updated LLR
		)
	{
	Run(tics,(const double*)r,(int*)llri,(int*)llro);
	return; 
	}

#endif

	//! Size of frame with pilot
	int K;
	
	//! Data frame size
	int Kdata;			

private:
	double afR,afI;
	const Modulator* mod;
	double* pil;
	double* buff_ab;
	double *absc;
	double ene;
	double *c;
	double sigmaD;
	double sigma;
	double fact;
	inline double LogBessSqrt(double x);
};

