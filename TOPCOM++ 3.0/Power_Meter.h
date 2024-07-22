// Power_Meter.h: interface for the Power_Meter class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdio.h>
#include "Include.h"

/*! \file 
\brief Declaration of the class Power_Meter
*/
/*! \ingroup Measurements 
\brief  Measures the power of a digital sequence.


The Power_Meter class measures the
average power \f$P_s\f$ of the input signal 
\f$s_{in}[n]\f$:
\f[
P_s=\frac{1}{N}\sum_{n=1}^{N}|s_{in}[n]|^2.
\f]
Through the method SetParameters() the user can specify the type of
input signal (real or complex). The Reset() method is used to reset
the power meter end the method Display() prints the results in an
output file defined by the user.

\author Guido Montorsi
*/

class Power_Meter  
{
public:
	Power_Meter();
	virtual ~Power_Meter();
	//! Set the main parameters of the block
	void SetParameters(const bool iscmplx=false	//!< Is the input signal complex? (default:false)
		);
	//! Run the power meter
	void Run(int ntics,			//!< Number of samples to be measured
		const double* input		//!< Signal to be measured
		);

	void Reset()
	{energy=0.;nsamples=0LL;};

	//! Get the current power
	double GetPower();

	//! Display the results
	void Display(FILE* stream=stdout		//!< Name of the output file (default:stdout)
		);
#ifdef CTOPCOM
	void Run(int ntics, const cmplx* input)
	{
		iscmplx=true;
		Run( ntics, (const double*)input);
	}
#endif
private:
	bool iscmplx;		//!< Is the input signal complex? 
	double energy;		//!< Energy
	__int64 nsamples;	//!< Number of samples
};

