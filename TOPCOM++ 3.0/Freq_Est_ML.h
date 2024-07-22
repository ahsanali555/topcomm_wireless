// Freq_Est_ML.h: interface for the Freq_Est_ML class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "FFT.h"
/*! \file 
\brief Declaration of class Freq_Est_ML
*/
/*! \ingroup Synchronization
\brief Maximum Likelihood Frequency Estimator (Rife and Boorstyn).

The main blocks of the algorithms are a block that perform the Fast
Fourier Transform of a set of input samples and an interpolating
algorithm to find the maximum of the absolute value of the transform.
The main parameters, initialized directly in the constructor, are the
base 2 logarithm of the size of input block $n$ and the ''pruning''
factor \f$n_2\f$.

A block of \f$N=2^n\f$ input samples fill the first part of a data block
of size \f$M=2^{n+n_2}\f$, the remaining part is zero padded.

FFT is performed on the block of \f$M\f$ complex data. The algorithm then
search for the maximum of the absolute value of the FFT and refine
the abscissa yielding the maximum, corresponding to the current
frequency estimate, by an interpolation using the trapezoidal rule.

Each tic of the Run() method process a single complex input. Input
are collected until a new block of \f$N\f$ data is ready to perform a new
estimation of the frequency.
  \see
	- DelayAndMultiply
	- Freq_Est_Kay
	- Freq_Est_Fitz
	- Freq_Est_LR
	- Freq_Est_Tretter
	- Freq_Est_YC

		For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".
*/
class Freq_Est_ML  
{
public:

	//! Construction and initialization of main parameters
	Freq_Est_ML(
		const int n,   //!< Log2 of the observation window
		const int n2=2 //!< Log2 of the pruning factor
		);

	//! Simulation method */
	/*! A tic processes a single input sample and updates the frequancy estimate stored in ::f. 
	The samples are collected until a new block of data is ready. Frequency estimate is then updated
	periodically. If provided, the sequence of estimates is also written in the output buffer */ 
	void Run(
		const int tics,			//!< Number of processed samples
		const double* input,	//!< Input complex sammples
		double*wc=0				//!< Estimated frequency samples
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,			//!< Number of tics
		const cmplx* input,				//!< Input complex samples
		cmplx*wc=0						//!< Estimated frequency samples
		)
	{
		Run(tics,(const double*)input,(double*)wc);
		return;
	}
#endif

	double f;					//!< Current frequency estimate

	virtual ~Freq_Est_ML();

private:
	double* buffer;
	FFT *fft;
	int n;
	int n2;
	int N,M;
	int now;
};

