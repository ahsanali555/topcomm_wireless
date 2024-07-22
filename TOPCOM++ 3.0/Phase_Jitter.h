// Phase_Jitter.h: interface for the Phase_Jitter class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "Gaussian_Process.h"

/*! \file
\brief Declaration of the class Phase_Jitter and all its interfaces
 */
/** \ingroup Channels
\brief Generates a phase jitter process and applies it to a signal

 The Phase_Jitter class embeds an instance of an object Gaussian_Process which generates  the Gaussian correlated  
process that describes the jitter. This object
must be instantiated and initialized outside the class and passed as parameter to the class itself.

SetParameters() specifies the process that
describe the phase jitter  \f$\theta_n\f$ in radians. The output complex
signal is
\f[z_n = x_n \exp(j\theta_n)\f]

  \see 
  - Phase_Jitter_from_Mask
  - DVBS2_Phase_Jitter
  - DVBRCS_Phase_Jitter
  - XBand_Phase_Jitter

For an example of its use see e.g. the test program "test_phase_jitters.cpp".

\author Guido Montorsi
*/
class Phase_Jitter  
{
public:
	Phase_Jitter();
	virtual ~Phase_Jitter();

	//! Specify the process that describe the phase jitter \f$\theta_n\f$ in radians
	void SetParameters(
			Gaussian_Process* freqproc	//!< Pointer to a Gaussian_Process
			);

	//! Apply the phase jitter to an input signal
	/*! The output buffer can coincide with the input buffer */
	void Run(int tics,		//!< Number of samples of the complex signal
		const double* inp,	//!< Input (complex) signal 
		double* out,			//!< Output (complex) signal
		double* ph=0			//!< Optional sequence of phase generated
		);

	//! Store in the output buffer the generated phase jitter complex process
	void Run(int tics,		//!< Number of samples of the complex signal
		double* out			//!< Output complex phase jitter
		);

#ifdef CTOPCOM
	void Run(int tics,		//!< Number of samples of the complex signal
		const cmplx* inp,	//!< Input (complex) signal 
		cmplx* out,			//!< Output (complex) signal
		double* ph=0			//!< Optional sequence of phase generated
		)
	{
		Run( tics,	(const double*) inp,(double*) out,ph);
	}
	void Run(int tics,		//!< Number of samples of the complex signal
		cmplx* out			//!< Output complex phase jitter
		)
	{
		Run( tics,(double*) out);
	}
#endif

//! The embedded Gaussian process
Gaussian_Process *process;
private:
	static const int nstep;
	double *buf;

//	Spectrum *spectrum;
	friend
Phase_Jitter* Phase_Jitter_from_Mask(
	const int nfreq,
	const double *freq,
	const double* slopes,
	const double gainDC,
	const double Fs
	);

friend
Phase_Jitter* DVBS2_Phase_Jitter(
								 const int Fs
								 );

friend
Phase_Jitter* DVBRCS_Phase_Jitter(const int Fs
								  );

friend
Phase_Jitter* XBand_Phase_Jitter(const int Fs 
								  );

friend
Phase_Jitter* CCSDS_Phase_Jitter(const int Fs 
								  );

};

/*! \ingroup Interface
\brief Return a phase jitter with a  PSD mask specified by user */
Phase_Jitter* Phase_Jitter_from_Mask(
	const int nfreq,		//!< Number of frequencies for mask specification
	const double *freq,		//!< Frequencies for mask specification
	const double* slopes,		//!< Starting slope at specified frequencies
	const double gainDC,	//!< Gain at DC (dB)
	const double Fs			//!< Sampling Frequency
	);

/*! \ingroup Interface
\brief Phase jitter of the DVB-S2 standard.*/
Phase_Jitter* DVBS2_Phase_Jitter(
								 const int Fs //!< sample rate (1,5,10,25,50 Mbaud)
								 );

/*! \ingroup Interface
\brief  Phase jitter of the DVB-RCS standard.*/
Phase_Jitter* DVBRCS_Phase_Jitter(const int Fs //!< sample rate (1= 128 kbaud , 2 = 1Mbaud)
								  );

/*! \ingroup Interface
\brief  X-band phase jitter.*/
Phase_Jitter* XBand_Phase_Jitter(const int Fs  //!< sample rate 
								  );

/*! \ingroup Interface
\brief  CCSDS phase jitter. 
Taken from CCSDS green book 413.0-G-1, pag 3-47

		-1e1  -25  30 db/decade
		-1e3  -85  14.3 db/decade
		-5e3  -95  0  db/decade
		-2e5  -95  21.2568 dB/decade
		-3e6  -120 0 dB/decade
		-1e8  -120
*/
Phase_Jitter* CCSDS_Phase_Jitter(const int Fs  //!< sample rate 
								  );
