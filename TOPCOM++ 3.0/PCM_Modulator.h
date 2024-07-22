// PCM_Modulator.h: interface for the PCM_Modulator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of class PCM_Modulator 
 */
/*! \ingroup Modems
\brief PCM modulator.

The PCM_Modulator class implements a phase modulator with residual carrier.
The transmittes signal can be modeled as:
\f[
x(t)=\sqrt{2P_T}\sin\{[2\pi f_ct+mP(t)d(t)\},
\f]
where \f$P_T\f$ is the total power, \f$m\f$ is the {\em modulation index} 
in radians (\f$0<m<\pi/2\f$), \f$\omega_c\f$ is the carrier center frequency,
\f$d(t)\f$ is the binary data sequence with symbol rate \f$R_s=1/T_s\f$ and \f$P(t)\f$ 
is a square of sinewave subcarrier
The possible data formats for the binary sequence \f$d(t)\f$ are:
- NRZ-L: ``1" is represented by level A and ``0" is represented by level B.
- NRZ-M (differential encoding):  ``1" is represented by a change in level and ``0" is represented by a no change in level.
- SPL: ``1" is represented by level A during the first half-symbol followed by level B during the second half-symbol
 and ``0" is represented by level B during the first half-symbol followed by level A during the second half-symbol.


The user specifies:
- The number of samples per signal.
- The modulation index \f$m\f$
- The line code type (NRZL, NRZM, or SPL).
- The subcarrier waveform type (no subcarrier, square subcarrier,sine subcarrier).
- The subcarrier frequency normalized with respect to symbol rate \f$f_c/R_s\f$.
- The Subcarrier phase \f$\theta_{sc}\f$ (in radians), only for sine subcarrier.


For an example of its use see e.g. the test program "test_PCM.cpp".

\author Gabriella Bosco

*/

class PCM_Modulator  
{
public:
	PCM_Modulator();
	virtual ~PCM_Modulator();
	//! Set the main parametrs of the block
	void SetParameters(const int nsamples,					//!< Number of samples per signal
						const double modindex,				//!< Modulation index	
						const char linecodtype[]="NRZL",	//!< Line code type ("NRZL","NRZM","SPL") (default: NRZL)	
						const int sc_waveform=1,			//!< Subcarrier waveform (0=no subcarrier,1=square subcarrier,2=sine subcarrier) (default: 1)	
						const int sc_frequency=1,			//!< Subcarrier frequency (normalized with respect to symbol rate) (default: 1)	
						const double sc_phase=0.			//!< Subcarrier phase (in radians) (default: 0.)
						);
	//! Run the modulator
	void Run(const int ntics,				//!< Number of input bits
				const int* inp,				//!< Input bits
				double* out					//!< Output signal
	);	

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,				//!< Number of input bits
				const int* inp,				//!< Input bits
				cmplx* out					//!< Complex Output signal
	)
	{
	Run(ntics,inp,(double*) out);
	return;
	}
#endif

	friend class PCM_Demodulator;

private:
	//! Select the subcarrier type
	void SetSubCarrier(	const int sc_waveform,			//!< Subcarrier waveform (0=no subcarrier,1=square subcarrier,2=sine subcarrier)	
						const int sc_frequency=1,			//!< Subcarrier frequency (normalized with respect to symbol rate) (default: 1)
						const double sc_phase=0.			//!< Subcarrier phase (in radians) (default: 0)
							);
	//! Select the line code type
	void SetLineCode(const char linecodtype[]				//!< Line code type ("NRZL","NRZM","SPL")
					);
	bool diff;			//!< Flag for differential encoding
	int ns;				//!< Number of samples per bit
	double m;			//!< Modulation index		
	int buffer;			//!< Buffer for differential encoding

	double *cosp;		//!< Cosine of subcarrier
	double *sinp;		//!< Sine of subcarrier
	double *signal;		//!< Linecode

};

