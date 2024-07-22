// Freq_Est_MM.h: interface for the Freq_Est_MM class.
//
//////////////////////////////////////////////////////////////////////
#include "ctopcom.h"
#pragma once
/*! \file 
\brief Declaration of class Freq_Est_MM
*/
/*! \ingroup Synchronization
\brief Frequency Estimator based on Mengali & Morelli algorithm.

The estimator uses \f$ N\f$ estimates of the autocorrelation function
\f[
\hat{R}(m)=\frac{1}{L_0-m}\sum_{k=m}^{L_0-1}r_k r_{k-m}^* \;\;1\leq
m\leq N
\f]

To compute the estimated frequency offset as
\f[
\hat{f}_d =\frac{1}{2\pi T}\sum_{m=1}^N
w(m)\times[\arg\{\hat{R}(m)\}-\arg\{\hat{R}(m-1)\}]
\f]

where the weight \f$w(m)\f$ are given by:
\f[
w(m)=\frac{3(L_0-m)(L_0-m+1)-N(L_0-N)}{N(4N^2-6NL_0+3L_0^2-1)}
\f]

The parameters supplied to the constructor are the  length of the
observation window \f$L_0\f$, and the number of estimated samples of
autocorrelation function \f$N\f$.

Each tic of the Run() method processes a complex input sample and
updates the frequency estimate ::f, a public variable of the class,
which is also optionally written in an output buffer specified by the
user.


  \see
	- DelayAndMultiply
	- Freq_Est_Fitz
	- Freq_Est_Kay
	- Freq_Est_LR
	- Freq_Est_ML
	- Freq_Est_Tretter
	- Freq_Est_YC

	For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".

*/
class Freq_Est_MM  
{
public:
	Freq_Est_MM(const int L0, const int N);
	virtual ~Freq_Est_MM();

	//! Simulation method. 
	/*! A tic processes a single input sample and updates the frequency estimate 
	stored in ::f. 
	if provided, the estimate is also written in the output buffer */ 
	void Run(const int tics,			//!< Number of tics
		const double* input,	//!< Input complex samples
		double *fest=0			//!< Sequence of frequency estimates [optional]
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,			//!< Number of tics
		const cmplx* input,				//!< Input complex samples
		cmplx*	fest=0					//!< Sequence of frequency estimates [optional]
		)
		{
		Run(tics,(const double*)input,(double*)fest);
		return;
	}

#endif

	double f;					//!< Current frequency estimate;

private:
	int L0;
	int N;
	double* R;
	double* line;
	double* w;
	int now;

};
