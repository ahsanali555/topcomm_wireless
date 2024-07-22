// OQPSK_Demodulator.h: interface for the OQPSK_Demodulator class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "Modulator.h"
#include "Demodulator.h"
#include "Sampler.h"
#include "OQPSK_Modulator.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class OQPSK_Demodulator 
*/
/*! \ingroup Modems
\brief OQPSK dedulator.
\author Gabriella Bosco

The OQPSK_Demodulator class is composed by  a delay line which eliminates the offset 
between the two quadratures,
followed by a FIR filter (defined as a member of teh Filter class), 
a sampler and a standard QPSK demodulator. Both hard or soft demodulation is available.

For an example of its use see e.g. the test program "test_OQPSK.cpp".
*/

class OQPSK_Demodulator  
{
public:
	OQPSK_Demodulator();
	virtual ~OQPSK_Demodulator();
	//! Set the main parameters of the block
	void SetParameters(const OQPSK_Modulator *mod,	//<! Reference Modulator
						const int iopt=0);		//<! Optimum sampling instant

	//! Set the variance of noise (for soft demodulation only)
	void SetSigma(
		const double sigma2,			//!< Noise variance
		const double fact=8.			//!< Noise variance
		);	
	
	//! Run the demodulator
	void Run(const int ntics,					//!< Number of symbols
		const double* inp,						//!< Input complex signal
		int* out,								//!< Output bits
		bool soft=false							//<! If soft=true, soft demodulation is used
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
		void Run(const int ntics,					//!< Number of symbols
		const cmplx* inp,						//!< Complex Input signal
		int* out,								//!< Output bits
		bool soft=false							//<! If soft=true, soft demodulation is used
		)
		{
		Run(ntics,(const double*) inp,out,soft);
		return;
		}
#endif

private:
	Demodulator *Demod;			//<! QPSK Demodulator
	Filter *Filt;				//<! Receiver filter
	Delay *Offset;				//<! Delay block to eliminate the offset between the two quadratures
	Sampler *Samp;				//<! Sampler
	double *sampled;			//<! Sampled signal
	double *filtered;			//<! Filterred signal
	int noffset;				//!< Offset (one half of the number of samples per signal) 
	int ns;						//!< Number of samples per signal 
	double *aligned;			//!< Received signal with the two quadratures aligned
};

