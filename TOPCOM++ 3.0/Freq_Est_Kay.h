// Freq_Est_Kay.h: interface for the Freq_Est_Kay class.
//
//////////////////////////////////////////////////////////////////////
#include "ctopcom.h"
#pragma once
/*! \file 
\brief Declaration of class Freq_Est_Kay
*/
/*! \ingroup Synchronization
\brief Frequency Estimator based on Kay algorithm.

The algorithm is basically a filtered delay and multiply with a
linear weighting function
\f[
\hat{f}_k=\sum_{i=1}^{N-2} w_i \Delta_{k-i}
\f]

where \f$\Delta_k=\tan^{-1}\left(\frac{x_{k,I}}{x_{k,R}}\right)\f$, \f$c_k\f$
are the received samples, and the coefficients are computed as

\f[
w_i =
\frac{3/2N}{N^2-1}\left(1-\left[\frac{t-(N/2-1)}{N/2}\right]^2\right)
\f]

Each tic of the Run() method process a single complex input and
update the frequency estimate f, which is also optionally written
in an output buffer.

  \see
	- DelayAndMultiply
	- Freq_Est_ML
	- Freq_Est_Fitz
	- Freq_Est_LR
	- Freq_Est_Tretter
	- Freq_Est_YC

		For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".


*/
class Freq_Est_Kay  
{
public:
	//! Construction and main parameter initialization*/
	Freq_Est_Kay(
		const int N		//!< Number of observed samples
		);

	//! Simulation method. 
	/*! A tic processes a single input sample and updates the frequancy estimate stored in ::f. 
	if provided, the estimate is also written in the output buffer */ 
	void Run(const int tics,			//!< Number of tics
		const double* input,		//!< Input complex samples
		double*fest=0				//!< Sequence of frequency estimates [optional]
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,			//!< Number of tics
		const cmplx* input,				//!< Input complex samples
		cmplx*fest=0					//!< Sequence of frequency estimates [optional]
		)
	{		
		Run(tics,(const double*)input,(double*)fest);
		return;
	}
#endif

	double f;					//!< Current frequency estimate;

	virtual ~Freq_Est_Kay();

private:
	double* w;
	double* line;
	int N;
	int now;
	double old0,old1;
};
