#pragma once
#include "DVBRCS_Transmitter.h"
#include "Demodulator.h"
#include "Interleaver.h"
#include "PN_Source.h"
#include "Filter.h"
#include "Delay.h"
#include "SISO_Decoder.h"
#include "CPM.h"
#include "Puncturer.h"
#include "Sampler.h"
#include "ESA/TurboPHI_Decoder.h"
#include "Spectrum.h"
#include "ShiftFrequency.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif


/*! \file
\brief Declaration of the class DVBRCS_Receiver
*/

/*! \ingroup Systems

\brief DVB-RCS receiver system as specified in standard ETSI EN 301 545-2 V1.2.0 (2013-12)

The RX system includes the following blocks:

CCCPM scheme:
- [optional] Filtering to limit ACI
- Iterative decoding between CPM SISO detector and SISO convolutional decoder
- [optional] Pilot removal

TCLM scheme:
- [optional] Shaping with SRRC filter of variable roll-off
- [optional] Pilot removal
- Demodulation
- Deinterleving
- Iterative decoding

All optional blocks are inserted depending on the correspondent configuration of the DVBRCS_Transmitter, which
is passed to the block through the method SetParameters().
The transmitter must be then configured before the construction of the receiver.

\sa DVBRCS_Transmitter

\author Gabriella Bosco
*/

class DVBRCS_Receiver
{
public:
	//! Empty constructor
	DVBRCS_Receiver(void);
	//! Empty destructor
	~DVBRCS_Receiver(void);

	//! Set the pointer to the reference transmitter and the number of decoder iterations
	/*! 
	The pointer is used to get all the relevant transmitter parameters to properly construct
	and configure the corresponding receiver.
	The method has to called each time a modification of the transmitter occurs, in particular 
	when changing the MODCOD.
	*/
	void SetParameters(
		const DVBRCS_Transmitter* e,    //!< Reference transmitter 
		const int niterin				//!< Number of iterations between CPM SISO detector and SISO convolutional decoder
		);
	
	//! Turn on or off the stopping rule in the iterative decoder (for CCCPM)
	void SetStop(const bool stopin=true){stop=stopin;}

	//! Insert the SRRC shaping filter (for TCLM)
	int AddRXFilter(
		const int ns,				//!< Number of samples per symbol (ns=1 means no filter)
		const double rolloff,		//!< Roll-off of the filter
		const int Nfil=20			//!< Number of filter taps
		);

	//! Insert the Anti-ACI filter (for CCCPM)
	int AddAAFilter(
		const int ns,				//!< [in] Number of samples per symbol (ns=1 means no filter)
		const double rolloff,		//!< [in] Roll-off of the filter
		const double BWnorm=1.,		//!< [in] Normalized bandwidth (w.r.t. symbol time)
		const int Nfil=20			//!< [in] Number of filter taps
		);

	//! Set the main RX parameters
	void TuneParameters_CCCPM(int nfilters_in,     //!< Number of filters in the SISO front-end
								int nini_in,						//!< Number of trellis steps for SISO initialization
								int ngroup_in						//!< Number of trellis steps for SISO update
								);

	//! Run the DVB-RCS Receiver
	void Run(
		const int tics,				//!< [in] Number of processed frames
		const double* input,		//!< [in] Input samples (complex)  
		int* dec					//!< [out] Output decoded bits or soft bits
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	//! Run the DVB-RCS Receiver
	void Run(
		const int tics,				//!< [in] Number of processed frames
		const cmplx* input,			//!< [in] Input samples (complex)   
		int* dec					//!< [out] Output decoded bits or soft bits
		)
	{
	Run(tics,(double*) input,dec);
	return;
	}
	
#endif

	//! Set the value of the variance of noise to be used in the computation of the LLRs
	void SetNoise(double snr);

private:
	
	CPM_Demodulator *CPMDem;//!< Embedded  CPM demodulator. (only for CCCPM)
	SISO_Decoder* iSISO;	//!< Embedded inner SISO processor. (only for CCCPM)
	SISO_Decoder* oSISO;	//!< Embedded outer SISO processor. (only for CCCPM)
	Depuncturer* Depunct;	//!< Embedded Depuncturer (only for CCCPM)
	ShiftFrequency * shiftf; //!< Frequency shift for tilt re-insertion (only for CCCPM)

	Filter* RXFil;			//!< Embedded SRRC filter  (only for TCLM)
	Demodulator* LinDem;	//!< Embedded demodulator  (only for TCLM)
	TurboPHI_Decoder* Dec;  //!< Embedded TurboPHI decoder (only for TCLM)
	Sampler *Samp;		    //!< Embedded sampler (only for TCLM)
	Delay *Align;		    //!< Embedded delay (only for TCLM)
	Spectrum *PS;

	double fact;			//!< factor for quantization of LLR
	int niter;				//!< Number of iterations SISO decoder

	double	*buffd1;		//!< Internal temporary buffer for doubles
	double	*buffd2;		//!< Internal temporary buffer for doubles
	int 	*buffi1;		//!< Internal temporary buffer for integers
	int 	*buffi2;		//!< Internal temporary buffer for integers

	int *llr;		//!< Temporary llr
	int *ext1;		//!< Temporary llr
	int *ext2;		//!< Temporary llr
	int *ext3;		//!< Temporary llr
	int nti;		//!< Number of trellis steps inner
	int nto;		//!< Number of trellis steps outer
	
	int nsymbols;	//!< Number of symbols
	bool stop;

	void SetParameters_CCCPM(const DVBRCS_Transmitter* e,	//!< [in] Reference tranmsitter
								  const int nfilters,
								  const int nini,
								  const int ngroup,
								  const int niterin,		//!< [in] Number of iterations between CPM SISO detector and SISO convolutional decoder
								  const double factor);

	void SetParameters_TCLM(const DVBRCS_Transmitter* e,	//!< [in] Reference tranmsitter
								  const int niterin			//!< [in] Number of iterations between CPM SISO detector and SISO convolutional decoder
								  );		

	void Run_CCCPM(
		const int tics,			//!< [in] Number of processed frames
		const double* input,	//!< [in] Input observations 
		int* dec				//!< [out] Output decoded bits or soft bits
		);


	void Run_TCLM(
		const int tics,			//!< [in] Number of processed frames
		const double* input,	//!< [in] Input observations 
		int* dec				//!< [out] Output decoded bits or soft bits
		);



	const DVBRCS_Transmitter* TX;	//!< Pointer to the transmitter storing the configuration
	
	int ns;					//!< Number of input samples per symbol
	double rolloff;			//!< Roll-off factor of RC filter
	int Nfil;				//!< Number of taps of RC filter
	int nfilters;			//!< Number of filters in the SISO front-end
	int nini;				//!< Number of trellis steps for SISO initialization
	int ngroup;				//!< Number of trellis steps for SISO update

};
