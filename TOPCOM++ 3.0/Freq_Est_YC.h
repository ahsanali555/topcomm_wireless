// Freq_Est_YC.h: interface for the Freq_Est_YC class.
//
//////////////////////////////////////////////////////////////////////
#include "ctopcom.h"
#pragma once
// WIN32 > 1000

/*! \file 
\brief Declaration of class Freq_Est_YC
*/
/*! \ingroup Synchronization
\brief Frequency estimator based on the Young & Crozier algorithm

   The parameters supplied to the constructor are the exponent of the input nonlinearity, the
   length of the observation window, the number of branches and the sorted vector
   with delay of each branch.
\author Guido Montorsi
\see
	- DelayAndMultiply
	- Freq_Est_Kay
	- Freq_Est_ML
	- Freq_Est_Fitz
	- Freq_Est_LR
	- Freq_Est_YC

For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".

*/


class Freq_Est_YC  
{
public:
	Freq_Est_YC(
		const int M,	//!< Exponent of input non linearity
		const int K,	//!< Observation window (number of samples) 
		const int B,	//!< Number of branches
		const int* delay //!< Sorted vector with delay of each branch
		);
	virtual ~Freq_Est_YC();

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
		cmplx*fest=0					//!< Sequence of frequency estimates [optional]
		)
	{
		Run(tics,(const double*)input,(double*)fest);
		return;
	}
#endif

	double f;					//!< Current frequency estimate;
	double q;					//!< Current phase estimate;

private:
	int M;	//!< Exponent of input non liearity
	int K;	//!< Observation window (number of samples) 
	int B;	//!< Number of branches
	int* delay; //!< Sorted vector with delay of each branch
	int maxd;
	double* line;
	double* acc;
	int time;

};

