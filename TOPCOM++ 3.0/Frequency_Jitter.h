// Frequency_Jitter.h: interface for the Frequency_Jitter class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "Gaussian_Process.h"
/*! \file
\brief Declaration of the class Frequency_Jitter
*/
/*! \ingroup Channels
\brief Generates a frequency jitter process and apply it to a signal

  The Frequency_Jitter class
embeds an instance of an object Gaussian_Process which generates  
the Gaussian correlated  
process that describes the jitter. This object
must be instantiated and initialized outside the class and passed as
 parameter to the class itself.

  The method SetParameters() allows to specify the process that
describe the frequency jitter \f$f_n\f$. The output complex signal is
\f[ z_n =  x_n \exp(j\theta_n)\f]
where
\f[\theta_n=\sum_{k=0}^n f_k.\f]

For an example of its use see e.g. the test program "test_jitters.cpp".

\author Guido Montorsi
*/
class Frequency_Jitter  
{
public:
	Frequency_Jitter();
	virtual ~Frequency_Jitter();
	//! Specify the process that describe the frequency jitter \f$f_n\f$
	void SetParameters(
			Gaussian_Process* freqproc  //!< Pointer to a Gaussian_Process
			);
	//! Applies the frequency jitter to an input signal
	/*! The output buffer can coincide with the input buffer */
	void Run(int tics,		//!< Number of samples of the input complex signal
		const double* inp,	//!< Input (complex) signal 
		double* out			//!< Output (complex) signal
		);
#ifdef CTOPCOM
	void Run(int tics,		//!< Number of samples of the input complex signal
		const cmplx* inp,	//!< Input (complex) signal 
		cmplx* out			//!< Output (complex) signal
		)
	{
		Run(tics,		
			(const double*) inp,	
			(double*) out			
			);
	}
#endif
private:
	Gaussian_Process *process;
	static const int nstep;
	double *buf;
	double phase;
};


