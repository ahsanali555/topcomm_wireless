// Envelope.h: interface for the Envelope class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdio.h>
#include <stdlib.h>

/*! \file 
\brief Declaration of the class Envelope
*/
/*! \ingroup Measurements 
\brief  Compute the absolute value of a complex signal.

For an example of its use see e.g. the test program "test_measurements.cpp".

\author Gabriella Bosco
*/

class Envelope  
{
public:
	Envelope(){};
	virtual ~Envelope(){};
	//! Run the block
	void Run(const int ntics,		//!< Number of signal samples
			const double* Input,	//!< Input signal
			double* Output			//!< Output signal
		);
#ifdef CTOPCOM
	void Run(const int ntics,		//!< Number of signal samples
		const cmplx* Input,	//!< Input signal
		cmplx* Output			//!< Output signal
		)
	{
		Run( ntics,(const double*) Input,(double*) Output);
	}
#endif
private:

};

