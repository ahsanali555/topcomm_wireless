#pragma once
#include "DVBSX_Transmitter.h"
#include "Demodulator.h"
#include "LDPC_Decoder.h"
#include "BCH_Decoder.h"
#include "Interleaver.h"
#include "Filter.h"
#include "Delay.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of the class DVBSX_Receiver
*/

/*! \ingroup Systems

\brief DVB-S receiver system as specified in standard ETSI EN 302 307 V1.3.1 (2013-03)
with modification introduced in DVB-SX

The RX system includes the following blocks:
- [optional] Shaping with SRRC filter of variable roll-off and decimation
- [optional] Fractionally spaced equalizer with ACG and phase estimation (see SetDetector())
- [optional] Physical Layer unscrambling
- [optional] Pilot stripping
- Computation of LLR on bits 
- [optional] Deinterleaving
- [optional] LDPC Iterative decoder
- [optional] BCH Decoder 
- [optional] Descrambler

all optional blocks are inserted depending on the correspondent configuration of the DVBSX_Transmitter, which
is passed to the block through the method SetTransmitter().
The transmitter must be then configured before the construction of the receiver.

\sa DVBSX_Transmitter

\author Guido Montorsi
*/
class DVBSX_Receiver
{
public:

	DVBSX_Receiver(void);
	~DVBSX_Receiver(void);

	//! Set the pointer to the reference DVB-SX transmitter
	/*! 
	The pointer is used to get all the relevant transmitter parameters to properly construct
	and configure the corresponding receiver.
	The method has to called each time a modification of the transmitter occurs, in particular 
	when changing the MODCOD.
	*/
	int SetParameters(
		const DVBSX_Transmitter*);

	//! Optionally change the precision for the fixed point LLR representation (default value set to 8)
	void SetLLRPrecision(
		const double fact //!< factor used before truncation: LLRq= NINT[LLR*fact] 
		);
	
	//! Activate the complete detector 
	/*! 
	The complete detector includes 
	- Adaptive Fractionally spaced equalizer
	- Phase compensation algorithm
	The usage of the complete detector is possible only if pilot, frame header and shaping filter
	are inserted in the DVBSX_Transmitter.
	*/
	void SetDetector(const double aeq,const int Nf);	

	//! Returns a boolean if the complete detector is trained.
	/*! This state is determined based on the values of the updating step of the fractionally spaced equalizer,
	the number of taps of the equalizer and the current number of performed updating steps. 
	*/ 
	bool IsTrained()
	{
		if(nup*aeq>2.*Nf*nsf)return true;
		else		 return false;
	}

	//! Set the relative delay of pilot sequence w.r.t. the observed samples 
	void SetDelay(
		const int delay //!< number of samples of delays
		);

	//! Displays the current values of Fractionally Spaced equalizer taps.
	/*! 
	 Nothing is diplayed if the detector block is not active.
	*/
	void DisplayFStaps(FILE* file=stdout);

	//! Reset centroid statistic
	/*!
	When the user passes the pointer to the indexes of the transmitted points through the optional
	parameter \ref in the method Run(), the detector carries on a running statistic of centroids.
	The running average  and variance of each centroid are then available and can be used to estimate the 
	effective signal to noise ratio.
	The centroids can also be used as internal reference constellation points instead of the nominal one
	(SetCentroids()). The running averages are computed using a one pole filter with tunable time constant.
	*/
	void ResetCentroids(double ace=1e-4 //!< Time constant for the centroids running averages
		);	

	//! Set the current centroid statistic  as demodulator reference 
	/*! 
	See ResetCentroids()
	*/
	void SetCentroids();						

	//! Displays the centroids and their variances  
	void DisplayCentroids(FILE*file=stdout);

	//! Displays the estimate of SNR based on centroid statistic
	/*! The SNR is evaluated as the ratio between the average energy of centroids and the average variance. 
	*/
	void DisplayCentroidSNR(FILE*file=stdout);	

	//! Run the DVB-S2/X Receiver
	int Run(
		const int frames,	//!< [in] Number of processed DVB-SX frames
		double *rx,			//!< [in] Input samples (complex) 
		int *dataRX,		//!< [out] Output decoded bits or soft bits
		const int *data=0,	//!< [in] Optional pointer to the transmitted bits for stopping LDPC decoder iterations
		const int *ref=0	//!< [in] Optional pointer to the indexes of transmitted points
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(
		const int frames,	//!< [in] Number of processed DVB-SX frames
		cmplx *rx,			//!< [in] Input samples (complex) 
		int *dataRX,		//!< [out] Output decoded bits or soft bits
		const int *data=0,	//!< [in] Optional pointer to the transmitted bits for stopping LDPC decoder iterations
		const int *ref=0	//!< [in] Optional pointer to the indexes of transmitted points
		)

	{
	return(Run(frames,(double*)rx,dataRX,data,ref));
	}
	
#endif

	int optdec;		    //!< Flag for optimal (0), suboptimal (min-sum-offset) (1), or layered suboptiomal (2)  LDPC decoding	
	bool dataequal;		//!< Flag to activate the usage of data for tuning the FS equalization
	int niter;			//!< Maximum number of iterations of LDPC decoder
	Demodulator* Dem;	//!< Embedded demodulator	
	FILE* filescat;		//!< Pointer to a file to store the scattering diagram


	double Pav;				//!< Current estimation of power at the output of fractionally spaced equalizer
	double MSE;				//!< Current estimation of MSE at the output of fractionally spaced equalizer
	double aeq;				//!< Time constant for updating Fractionally Spaced  equalizer and AGC

	LDPC_Decoder*	Dec;		//!< Embedded LDPC decoder	
	BCH_Decoder*	BCHDec;		//!< Embedded BCH decoder
private:

	Modulator* LocMod;			//!< Embedded local modulator with used constellation points	
	Filter* RXFil;			//!< Embedded SRRC filter
	double fact;			//!< factor for quantization of LLR
	int time;		//!< Internal clock
	int	*buffi1;		//!< Internal temporary buffer for integers
	int	*buffi2;		//!< Internal temporary buffer for integers
	double	*buff1;		//!< Internal temporary buffer for doubles
	double	*buff2;		//!< Internal temporary buffer for doubles
	double* buffer;

	double* bufth;	//< Buffer for the phase estimate
	int timeth;		//< pointer to the phase buffer
	int gap;		//< Gap from last phase estimate



	const DVBSX_Transmitter* TX;	//!< Pointer to the transmitter storing the configuration
	//! Physical layer scrambler
	int xi;			//! < sco
	int yi;
	int nsc;		//!< Number of Gold sequence for PL scambler
	int ns;
	int rate;
	int blen;
	double rolloff;
	int Nfil;

// Storage for complete detector

	double RunDetector(const double* inp, double* out, const int* ref=0);
	double theta;			//!< Current phase estimate
	double* f;				//!< Equalizer coefficients
	double* linef;			//!< Equalizer line
	double* pil;		//!< Reference pilot sequence

	int Nsymb;				//!< Periodicity of pilot sequence
	Delay* del;			//!< Delays to align sequence to output of phase estimation
	int Nf;			//!< Number of filter taps
	int nowf;		//!< Pointer to input of equalizer
	int EQlength;	//!< Length of FS equalizer line
	int nsf;		//!< sample per symbol FS equalizer
	double acct[2];	//!< Accumulators of real and imaginary parts for phase estimation on pilots
	double thn;		//!< Accumulator for phase estimation on pilots
	double tho;		//!< Accumulator for phase estimation on pilots

	double* centroid;
	double* varcen;
	double ace;
	int nup;
	int beta;

	double PhaseInterp(const int index, const double* z);

	Delay* delref;
	Delay* delref2;
};
