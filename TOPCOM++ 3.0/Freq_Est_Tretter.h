// Freq_Est_Tretter.h: interface for the Freq_Est_Tretter class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "ctopcom.h"


/*! \file 
\brief Declaration of class Freq_Est_Tretter
*/
/*! \ingroup Synchronization
\brief Frequency estimator based on the Tretter algorithm

   The only parameter supplied to the class with the constructor 
   is the number of observation samples.
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

class Freq_Est_Tretter  
{
public:
	//! Construction and main parameter initialization*/
	Freq_Est_Tretter(
		const int N		//!< Number of observed samples
		);

	//! Simulation method. 
	/*! A tic processes a single input sample and updates the frequancy estimate 
	stored in ::f. 
	if provided, the estimate is also written in the output buffer */ 
	void Run(const int tics,			//!< Number of tics
		const double* input,	//!< Input complex samples
		double*fest=0			//!< Sequence of frequency estimates [optional]
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,			//!< Number of tics
		const cmplx* input,				//!< Input complex samples
		cmplx* fest=0					//!< Sequence of frequency estimates [optional]
		)
	{
		Run(tics,(const double*)input,(double*)fest);
		return;
	}
#endif

	double f;					//!< Current frequency estimate;
	double q;					//!< Current phase estimate;

	virtual ~Freq_Est_Tretter();
private:
	int N,M;
	double K1;
	double*line;
	double oldy,a,b,c;
	int now;
};
