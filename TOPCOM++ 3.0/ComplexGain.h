// ComplexGain.h: interface for the ComplexGain class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <math.h>
/*! \file 
\brief Declaration of the class ComplexGain
*/
/*! \ingroup DSP 
\author Gabriella Bosco
\brief  Multiplies a complex signal by a complex constant.

For an example of its use see e.g. the test program "test_miscellanea.cpp".
*/
class ComplexGain  
{
public:
	//! Constructor: the complex gain is expressed in terms of real and imaginary part
	ComplexGain(
		const double* gain				   //!< Complex gain (real and imaginary part)
	);
	//! Constructor: the complex gain is expressed in terms of modulus and phase
	ComplexGain(
		const double mod,	//!< Modulus of the complex gain
		const double phase  //!< Phase of the complex gain (in radians)
	);
	//! Destructor
	virtual ~ComplexGain();
	
	//! Multiplies the input complex signal for a complex constant
	/*! The output buffer can coincide with the input buffer */
	void Run(const int ntics,		//!< Number of input signal samples
			const double* Input,	//!< Input signal
			double* Output			//!< Output signal
			);

	#if CTOPCOM
	void Run(const int ntics,		//!< Number of input signal samples
			const cmplx* Input,	//!< Input signal
			cmplx* Output			//!< Output signal
			)
	{
		Run(ntics,		//!< Number of input signal samples
				(const double*) Input,	//!< Input signal
				(double*) Output			//!< Output signal
				);
		return;
	}
	ComplexGain(const cmplx* gain)
	{
		ComplexGain((const double*)gain);
	}
	#endif

	void SetPolar(const double mod,		//!< Number of input signal samples
			const double phase	//!< Input signal
			){re=mod*cos(phase);im=mod*sin(phase);}

	double re,		//!< Real part of complex gain
		im;			//!< Imaginary part of complex gain
private:
};
