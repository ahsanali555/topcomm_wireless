// ComplexMultiplier.h: interface for the ComplexMultiplier class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class ComplexMultiplier
*/
/*! \ingroup DSP 
\brief  Multiplies two complex signals.

For an example of its use see e.g. the test program "test_miscellanea.cpp".

*/

class ComplexMultiplier  
{
public:
	ComplexMultiplier();

	virtual ~ComplexMultiplier();
	//! Multiplies two complex signals
	/*! The output buffer can coincide with one of the input buffers */
	void Run(const int ntics,		//!< Number of signal samples
			const double* Input1,	//!< First input signal
			const double* Input2,	//!< Second input signal
			double* Output		//!< Output signal
			);

	//! Multiplies a complex signals with a phase signal
		void RunPhase(const int ntics,		//!< Number of signal samples
			const double* Input1,	//!< First input signal
			const double* Phase,	//!< Phase rotation
			double* Output		//!< Output signal
			);

		
#if CTOPCOM
		void Run(const int ntics,		//!< Number of signal samples
			const cmplx* Input1,	//!< First input signal
			const cmplx* Input2,	//!< Second input signal
			cmplx* Output		//!< Output signal
			)
		{
			Run(ntics,(const double*) Input1,(const double*) Input2,(double*) Output);
			return;
		}

	//! Multiplies a complex signals with a phase signal
		void RunPhase(const int ntics,		//!< Number of signal samples
			const cmplx* Input1,	//!< First input signal
			const cmplx* Phase,	//!< Phase rotation
			cmplx* Output		//!< Output signal
			)
		{
			Run(ntics,(const double*) Input1,(const double*) Phase,(double*) Output);
			return;
		}
#endif
};
