// PCM_Demodulator.h: interface for the PCM_Demodulator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "PCM_Modulator.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class PCM_Demodulator 
*/
/*! \ingroup Modems
\brief PCM demodulator.

The PCM_Demodulator class implements a maximum-likelihood phase demodulator (both hard ans soft).
If the line code type is ``NRZ-L", differential decoding is performed.
The user specifies the variance of AWGN noise (used in soft decoding) and a factor setting 
the precision of quantization in soft decoding (see the Demodulator class for explanation).

For an example of its use see e.g. the test program "test_PCM.cpp".

\author Gabriella Bosco
*/
class PCM_Demodulator  
{
public:
	PCM_Demodulator();
	virtual ~PCM_Demodulator();
	//! Set the main parameters of the block
	void SetParameters(const PCM_Modulator* mod			//!< PCM modulator
		);
	//! Set the variance of noise for soft decoding
	void SetSigma(const double sigma2,		//!< Variance of AWGN noise
				  const double factin=8.	//!< Factor setting the precision of quantization in soft decoding.
					);
	//! Set the precision of llrs representation
	void SetPrec(const double fact		//!< Factor setting the precision of quantization in soft decoding.
					);
	//! Run the demodulator
	void Run(const int ntics,	//! Number of bits
		const double* inp,		//! Input signal
		int* out,			//! Output bits
		bool soft=false	//!< If soft=true, soft decoding is performed.
			);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
		void Run(const int ntics,	//! Number of bits
		const cmplx* inp,		//! Complex Input signal
		int* out,			//! Output bits
		bool soft=false	//!< If soft=true, soft decoding is performed.
			)
		{
		Run(ntics,(const double*) inp,out,soft);
		return;
		}
#endif


	private:
		//! Look-up table for max* operation
		inline int lut(int a, const int b);
		const PCM_Modulator *refmod;	//!< PCM modulator
		double fact;					//!< Conversion factor
		double a;						//!< a=fact/2/sigma2
		int *ic;						//!< Look-up table for max* operation
		int maxtab;						//!< Dimension of the look-up table
		int softbuf;					//!< Buffer used for differential soft decoding
		int buf;						//!< Buffer used for differential hard decoding
	
};
