// OQPSK_Modulator.h: interface for the OQPSK_Modulator class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "Modulator.h"
#include "Filter.h"
#include "Delay.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class OQPSK_Modulator 
*/
/*! \ingroup Modems
\brief OQPSK modulator.

The OQPSK_Modulator class implements a  modulator for offset-QPSK (OQPSK) signals.
An OQPSK modulation consists of a standard QPSK modulation in which the two quadrature components are
offset in time by a bit period \f$T_B\f$. 

The OQPSK modulator is composed by  a standard QPSK modulator, followed by a shaping filter
and a delay line which introduces the offset between the two quadratures. 
The user specifies:
- The number of samples per signal \f$n_s\f$ (which has to be equal to at least 10 in order 
to obtain reliable results). 
- The shaping filter, defined as  a member of the class Filter.

For an example of its use see e.g. the test program "test_OQPSK.cpp".

\author Gabriella Bosco
*/

class OQPSK_Modulator  
{
public:
	OQPSK_Modulator();
	virtual ~OQPSK_Modulator();

	//! Set the main parameters of the OQPSK modulator
	void SetParameters(const int nsamples,		//!< Number of samples per signal 
					   double alpha=0.5,		//!< Roll-off of SRRC shaping filter
					   int delay =-1		//!< Length of filter (number of symbols)
					   );

	//! Set the main parameters of the modulator
	void SetParameters(const int nsamples,		//!< Number of samples per signal 
					   Filter* filter=0		//!< Shaping filter
					   );

	//! Run the modulator
	void Run(const int ntics,		//!< Number of input symbols
				const int* inp,		//!< Input signal
					double* out		//! Output signal
					);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,	//!< Number of input symbols
		 const int* inp,		//!< Input signal
		 cmplx* out				//!< Complex Output signal
		 )
	{
	Run(ntics,inp,(double*) out);
	return;
	}
#endif

	friend class OQPSK_Demodulator;
	int delay;

private:
	Modulator *Modul;		//!< QPSK modulator
	Filter *Filt;			//!< Shaping filter
	Delay *Offset;			//!< Delay block to introduce the offset between the two quadratures
	double *mod;			//!< Modulator output
	int noffset;			//!< Offset (one half of the number of samples per signal) 
	int ns;					//!< Number of samples per signal 
	double *filtered;		//!< Filter output
};

