// Delay.h: interface for the Delay class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class Delay
*/
/*! \ingroup DSP 
\brief Insert a delay on a discrete sequence.

The block inserts an arbitrary delay over integer real and complex sequences.
Input buffer and output buffers can be the same

\author Guido Montorsi
*/

class Delay  
{
public:
	Delay();
	virtual ~Delay();

	//! Set the main parameter of the block
	void SetParameters(const int delay,				//!< Delay expressed in number of samples
		const int type				//!< Type of signal (0 - integer, 1 - real, 2 - complex double)
		);

	//! Run the block on integer signals 
	/*! The output buffer can coincide with the input buffer */
	void Run(const int tics,			//! Dimension of the input vector
		const int* input,			//! Input signal
		int* output					//! Ouput signal
		);

	//! Run the block on real or complex signals 
	/*! The output buffer can coincide with the input buffer */
	void Run(const int tics,			//! Number of input samples
		const double* input,	//! Input signal
		double* output,			//! Output signal
		const int diffdelay=0	//! If diffdelay=1 or 2, only real or imaginary part of input signal is delayed
		);

	//! Reset the delay buffer
	void Reset(const int init=0	//!< Initial values in delay
		);

#ifdef CTOPCOM
	void Run(const int tics,			//! Number of input samples
		const cmplx* input,	//! Input signal
		cmplx* output,			//! Output signal
		const int diffdelay=0	//! If diffdelay=1 or 2, only real or imaginary part of input signal is delayed
		)
	{
		Run( tics,			//! Number of input samples
			(const double*) input,	//! Input signal
			(double*) output,			//! Output signal
			diffdelay	//! If diffdelay=1 or 2, only real or imaginary part of input signal is delayed
			);
	}
#endif


private:
	int delay;			//!< Delay expressed in number of samples
	void* buffer;		//!< Buffer used to store the delayed samples
	int pointer;		//!< Time index
	int type;			//!< Type of signal (0 - integer, 1 - real, 2 - complex double)
};

