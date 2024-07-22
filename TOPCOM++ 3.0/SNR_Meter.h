// SNR_Meter.h: interface for the SNR_Meter class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdio.h>
#include "Include.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file
\brief Declaration of the class SNR_Meter
 */  

/** \ingroup Measurements

\brief  Measure the Signal to Noise Ratio (SNR) of a random variable

The SNR is defined as the ratio between the square of the mean and the variance of
the random variable.

Through the method Start_Stop() the user can specify the starting instant for the measure 
and the number of collected samples for the estimate.

The method Display() print the current SNR estimate on the specified output device.

\author Guido Montorsi
*/
class SNR_Meter  
{
public:
	SNR_Meter();
	virtual ~SNR_Meter();

	//! Set starting sample and number of samples considered
	void Start_Stop(
		const int istart=0,	 //!< Starting instant for measure. (from function call)
		const int nsmeas=-1	 //!< Number of samples considered. (-1: all)
		);

	//! Run the block
	void Run(
		const int tics,		//!< Number of samples.
		const double* inp	//!< Input samples
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
void Run(const int tics,   //!< Number of input samples
			const cmplx* inp		//!< Complex Input signal
			 )
	   {
	   Run(tics,(const double*) inp);
	   }	
#endif

	//! Display the results on the specidied output device.
	void Display(FILE* file=stdout);

	//! Reset the measure
	void Reset();

	double mean;	//!< Accumulated mean.
	double vqm;		//!< Accumulated mean square
	__int64 nsamp;	//!< Number of samples

private:


	bool start;
	int istart;
	int nsmeas;
};
