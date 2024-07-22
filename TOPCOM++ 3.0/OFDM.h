// OFDM.h: interface for the OFDM class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "FFT.h"
#include "Filter.h"
#include <stdio.h>
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declarations of OFDM modulator and demodulator */
/*! \ingroup Modems
\brief OFDM modulator.

 This block implements all the typical operations performed in a Orthogonal Frequency Division \Modulator.

  - Modulation of n carriers with arbitrary constellations
  - Complex IFFT
  - Addition of cyclic prefix samples.
  - Windowing.

  Additionally, the user can specify a pilot pattern and a channel allocation mechanism to
  divide the time frequency plane into logical channels and perform multiaccess emulation, see the manual
  for more details on these mechanisms.

For an example of its use see e.g. the test program "test_CompCodeComb.cpp".

  Interfaces: OFDM_DVBT()
\author Guido Montorsi
*/
class OFDM  
{
public:
	OFDM();

	virtual ~OFDM();

	//! Set the main parameters of OFDM modulator
	int SetParameters(const int n,//!< Log2 of the number of carriers
		const int Cp		  //!< Number of added samples for cyclic prefix.
		);

	//! Set the active (modulated) carriers 
	void SetActiveCarriers(
		const int* act	  //!< Vector of active carriers (1=active, 0=not active)
		);

	//! Set the pilot sequence 
	void SetPilots(
		const int  nsymb,		//!< Periodicity of pilot insertion pattern
		const int  per,			//!< Period of pilot sequence (number of OFDM symbols)
		const double* pilots,	//!< Pilot sequence
		const int npos,			//!< Number of pilot positions
		const int* pos			//!< Positions of pilots
		);		

	//! Display the frame structure 
	void Display(FILE* f=stdout);

	//! Set the logic channels associated to the modulator
	void SetChannels(
		const int perall,		//!< Allocation period (number of OFDM symbols)
		const int *alloc,		//!< Allocation vector
		const int gt=1,			//!< Granularity in time
		const int gf=1			//!< Granularity in frequency
		);

	//! Run the block with map of data for a single users
	int Run(const int nsymb,	//!< Number of generated OFDM symbols
		const double*input,	//!< Input complex  points  
		double* output,		//!< Modulated OFDM signal (with cyclic prefix)
		const int user=0	//!< Allocated user
		);

	//! Run the block with map of data of all users 
	int Run(const int nsymb,	 //!< Number of generated OFDM symbols
		const double**input, //!< Pointer to data of different users 
		double* output		 //!< Modulated OFDM signal (with cyclic prefix)
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
		//! Run the block with map of data for a single users
	int Run(const int nsymb,	//!< Number of generated OFDM symbols
		const cmplx* input,	//!< Input complex  points  
		cmplx* output,		//!< Modulated OFDM signal (with cyclic prefix)
		const int user=0	//!< Allocated user
		)
	{
	return (Run(nsymb,(const double*) input,(double*) output, user));
	}

	//! Run the block with map of data of all users 
	int Run(const int nsymb,	 //!< Number of generated OFDM symbols
		const cmplx**input, //!< Pointer to data of different users 
		cmplx* output		 //!< Modulated OFDM signal (with cyclic prefix)
		)
	{
		return (Run(nsymb,(const double**)input,(double*) output));
	}
#endif


	friend class OFDM_Demod;

	int n;	//!< Log-2 of number of carrers of OFDM
	int N;	//!< Number of carriers of OFDM
	int Cp;	//!< Length of cyclic prefix
	int Nac;		//!< Number of active carriers

private:
	double* buffer;
	FFT* fft;
	int* act;	//!<active carriers
	int* alloc;
	int* pos;	// Positions of pilot
	double* pilots;

	double  sq;
	int per;
	int npos;	// Pointer in position
	int nsymb;	// Frame length

	int time;//!< Internal clock
	int np1;		//Pointer in pilots pattern
	int np2;		//Pointer in pilots pattern

	// Channel allocation
	int perall;
	int gf;	//!< Frequency granularity
	int gt;	//!< Time granularity
	int N1;	 //!< Number of frequency slots
};

/**\ingroup Modems
\brief OFDM demodulator.

  The block is initialized from the correspondent OFDM modulator and performs all
  typical operations of an OFDM demodulator
  - Removal of cyclic prefix
  - FFT
  - Extraction of pilot symbols
  - [Optional Channel estimation]
 - Equalization 

 The output is the sequence of [equalized] data complex points, 
 normalized to have unitary average energy. These
 complex points can then be provided to the Demodulator for Hard or Soft decisions.

 For an example of its use see e.g. the test program "test_CompCodeComb.cpp".

 \author Guido Montorsi
*/
class OFDM_Demod 
{
public:
	OFDM_Demod();
	/** Constructor */
	OFDM_Demod(const OFDM*	//!< pointer to the correspondent modulator
	);

	virtual ~OFDM_Demod();

	/*! Set the main demodulator parameters */
	void SetParameters(
		const OFDM* mod //!< pointer to the correspondent modulator
		);


	//! Performs OFDM demodulation and optionally equalization 
	/*! On ouput, if required, it also return the squared magnitude of the channel gains for each sample. 
	This information should be used to proper compute the soft information in the following
	Demodulator
	*/
	int Run(
		const int nsamples,		//!< Number of input complex samples
		const double* input,	//!< Input complex samples 
		double* data,			//!< Demodulated data complex points (single user)
		double* chabsout=0,		//!< Optional buffer to store squared magnitude of channel gains
		int user=0				//!< Optional index of decoded user
		);

	/*! Activates channel estimation algorithm and equalization */
	void SetChannelEstimation(
		const double tau,		//!< Time constant for filtering in time
		const bool est  =true	//!< Activates (or deactivates) the channel estimation from pilots
		);

	/*! Print the current channel estimates */
	void DisplayChannelEstimation(FILE* file=stdout);

	/*! Pass to the demodulator the channel coefficients (frequency domain) */
	void SetChannelCoeff(
		const double *chcoeffin);


	/*! Current estimate of transfer function */
	bool estch;			//!< Flag that indicates that channel estimation is on
	bool equal;			//!< Flag that indicates if equalization is on
	double* chcoeff;	//!< Current estimated channel coefficients
	double* chabs;		//!< Current estimated amplitude of channel coefficients

private:
	const double *chideal;
	bool* chinterp;		//!<
	int np1,np2;
	FFT* fft;
	const OFDM *refmod;
	int time,no;
	double* buffer;//!< Internal buffer to perform FFT.

	// For channel estimation
	Filter** chfil;


friend
OFDM* OFDM_DVBT(const int mode,const int guard,const bool half);

friend
OFDM* OFDM_DVBT2(const int mode,const int PP,const int guard,const bool half);

friend
OFDM* OFDM_WiMAX_UL_PUSC(const int nfft, const int Cp, const int N_Tiles_used);
};



/*
*! \ingroup Interface
\brief OFDM modulator as specified in the DVBT standard.

This interface returns an OFDM modulator with features as specified
by the DVB-T standard. The only parameters are the mode
(0=2k or 1=8k) and the guard time (between 1 and 5) that specifies the
cyclic prefix. 1705 Carriers are active in the 2k mode and 6817 in
the 8k mode. Correct pilots are inserted according to the pattern
specified by the standard. As a result of the pilot overhead, the
number of available data slots per OFDM symbol is  1512 for the 2k
mode and 6048 in the 8k mode.
For half size modes (1k and 4k) set to true the parameter "half". 
\author Guido Montorsi
*/

OFDM* OFDM_DVBT(const int mode=0,	//!< Short (0=2k) or Long (1=8k) FFT 
				const int guard=5.,	//!< Specify cyclic prefix length (5=>32, 4=>16, ...,2=>4)
				const bool half=false //!< Half size of FFT (1k / 4k)
				);

/* \ingroup Interface
\internal
\brief OFDM modulator as specified in the DVB-T2 standard.

Work in progress.

\author Guido Montorsi
*/

OFDM* OFDM_DVBT2(const int mode,	//!< Log_2 of length FFT (1k=0,1,2,3,4,32k=5)  
				 const int PP,		//!< Pilot pattern (1..8)
				 const int guard=5.,	//!< Specify cyclic prefix length (5=>32, 4=>16, ...,2=>4)
				const bool half=false //!< Half size of FFT (1k / 4k)
				);


OFDM* OFDM_WiMAX_UL_PUSC(const int nfft, const int Cp, const int N_Tiles_used);