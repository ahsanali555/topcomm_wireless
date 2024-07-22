// Timing_Jitter.h: interface for the Timing_Jitter class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "AWGN_Channel.h"
#include "Gaussian_Process.h"
#include "AtoD.h"


/*! \file 
\brief Declaration of the class Timing_Jitter
*/
/*! \ingroup Channels
\author Guido Montorsi
\brief Generator of sampling jitter.

The Timing_Jitter class embeds an instance of an object Gaussian_Process which generates  the Gaussian correlated  
process that describes the jitter. This object
must be instantiated and initialized outside the class and passed as parameter to the class itself.
The generated Gaussian process, which is a zero
mean process, is offset by a constant and then fed to an AtoD object


The generated Gaussian process, which is a zero
mean process, is offset by a constant and then fed to an AtoD object 


SetParameters() specifies the two component blocks and optionally an average sampling frequency.
The average sampling frequency is added to generated zero-mean process.

The Timing jitter is generated applying the generated process \f$f_n\f$
as the instantaneous frequency of an AtoD
converter. The signal that describe the jitter is provided at the same sample rate of the
input signal \f$x_n\f$.

For an example of its use see e.g. the test program "test_jitters.cpp".

\author Guido Montorsi
*/
class Timing_Jitter  
{
public:
	Timing_Jitter();
	virtual ~Timing_Jitter();
	//! Specify the process that describes the istantaneous frequency of A/D
	void SetParameters(
		Gaussian_Process* freqproc,	//!< Timing jitter process
		AtoD* atodin,				//!< Analog to digital converter
		double avfreq=1.			//!< Average sampling frequency
		);

	//!  Specify the process that describes the istantaneous frequency of A/D.
	/*! Assumes a on pole filter for model of the istantaneous frequency process.*/
	void SetParameters(
		const double decay, //!< time constant
		const double pow	//!< power of process
		);

	//!  Specify the process that describes the istantaneous phase of A/D.
	/*! White A/D process */
	void SetParameters(
		const double stdev
		);

	//! Apply the Timing  jitter to an input signal
	int Run(int tics,	//!< Number of samples of the complex signal
		const double* inp,	//!< Input (complex) signal 
		double* out			//!< Output (complex) signal
		);
	AtoD* atod;			// Upsampler for timing jitter process

#ifdef CTOPCOM
	int Run(int tics,	//!< Number of samples of the complex signal
		const cmplx* inp,	//!< Input (complex) signal 
		cmplx* out			//!< Output (complex) signal
		)
	{
		Run(tics,(const double*) inp,(double*) out);
	}
#endif
private:

	Gaussian_Process *freqproc;
	double* frequency;
	double avfreq;
	double stdev;
	int seed;
	static const int nstep;
};