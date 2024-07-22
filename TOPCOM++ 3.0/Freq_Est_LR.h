// Freq_Est_LR.h: interface for the Freq_Est_LR class.
//
//////////////////////////////////////////////////////////////////////
#include "ctopcom.h"

#pragma once
/*! \file 
\brief Declaration of class Freq_Est_LR
*/
/*! \ingroup Synchronization
\brief Frequency Estimator based on Luise & Reggianini algorithm.

  User specifies the length of observation window and the number of branches 
  in the estimator.
  \see
	- DelayAndMultiply
	- Freq_Est_Kay
	- Freq_Est_ML
	- Freq_Est_Fitz
	- Freq_Est_Tretter
	- Freq_Est_YC

		For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".

*/

class Freq_Est_LR  
{
public:
	Freq_Est_LR(const int Nin, const int Min);
	virtual ~Freq_Est_LR();
	//! Simulation method. 
	/*! A tic processes a single input sample and updates the frequancy estimate stored in ::f. 
	If provided, the estimate is also written in the output buffer */ 
	void Run(const int tics,			//!< Number of tics
		const double* input,	//!< Input complex samples
		double*fest=0			//!< Sequence of frequency estimates [optional]
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

	double f;	//!< Current frequency estimate


private:
	int M;
	int N;
	double* line;
	double* w;
	double r0,r1;	//Sum of correlations
	int now;
};
