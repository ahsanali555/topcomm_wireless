// DelayAndMultiply.h: interface for the DelayAndMultiply class.
//
//////////////////////////////////////////////////////////////////////
#include "ctopcom.h"
#pragma once
/*! \file 
\brief Declaration of class DelayAndMultiply
*/
/*! \ingroup Synchronization
\brief Simple delay and multiply frequency estimator.

The frequency is estimated by taking the complex product of the current sample with the complex
conjugate of a sample in the past. 
The normalized phase of the obtained complex value is used as an estimator
of the frequency.
  
More sophisticated frequency estimator are
  \see
	- Freq_Est_Kay
	- Freq_Est_ML
	- Freq_Est_Fitz
	- Freq_Est_LR
	- Freq_Est_Tretter
	- Freq_Est_YC

	 For an example of its use see e.g. the test program "Frequency_Synchronization.cpp".

\author Guido Montorsi
*/
class DelayAndMultiply  
{
public:
	//! Construction and initialization
	DelayAndMultiply(
		const int delay	//!< Set the delays bewteen samples to multiply
		);
	virtual ~DelayAndMultiply();
	//! Simulation  method. 
	/*! One tics processes one input complex signal and produces one output 
	frequency estimate.	
	*/
	void Run(const int tics,	//!< Processed input samples
			const double* input,//!< Input samples
			double *f=0			//!< Sequence of frequency estimates [optional]
			);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,	//!< Processed input samples
			const cmplx* input, //!< Complex Input samples
			cmplx* f=0			//!< Sequence of frequency estimates [optional]
			)
	{
	Run(tics,(const double*) input,(double*)f);
	return;
	}

#endif

	double f;					//!< Current frequency estimate
private:
	double *line;
	int now;
	int   delay;
};
