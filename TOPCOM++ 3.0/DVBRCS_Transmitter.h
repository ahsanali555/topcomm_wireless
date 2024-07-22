#pragma once
#include "PN_Source.h"
#include "BCH_Encoder.h"
#include "LDPC_Encoder.h"
#include "ESA/TurboPHI_Encoder.h"
#include "Interleaver.h"
#include "Modulator.h"
#include "CPM.h"
#include "Filter.h"
#include "Transponder.h"
#include "Trellis_Encoder.h"
#include "Puncturer.h"
#include "ShiftFrequency.h"
	
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of the class DVBRCS_Transmitter
*/

/*! \ingroup Systems

\brief Transmitter system as specified in standard ETSI EN 301 545-2 V1.2.0 (2013-12)

The TX system includes the following blocks:

SCCCPM scheme:
- Convolutional encoding
- Puncturing
- Interleaving
- [optional] Pilot insertion
- CPM modulation

LMTC scheme:
- Turbo-phi encoding
- Interleaving
- [optional] Pilot insertion
- Linear modulation
- [optional] Shaping with SRRC filter of variable roll-off

All optional blocks can be inserted by using the methods for their configuration. 
The square root raised cosine transmitting filter is configured with the method AddTXFilter().
The insertion of known symbols is activated with the method AddKnownSymbols.

\sa DVBRCS_Receiver

\author Gabriella Bosco

*/
class DVBRCS_Transmitter
{
public:
	//! Empty constructor
	DVBRCS_Transmitter(void);
	//! Empty destructor
	~DVBRCS_Transmitter(void);
	
	//! Set the active MODCOD of DVB-RCS transmitter
	int SetMODCOD(
		const int CPM_flag,		//!< [in] Flag for CPM or linear modulation
		const int MODCOD,		//!< [in] Index of MODCOD
		const int ns			//!< [in] Number of samples per symbol
		);

	//! Insert the SRRC shaping filter
	int AddTXFilter(
		const int ns,			//!< [in] Number of samples per symbol (ns=1 means no filter)
		const double rolloff,	//!< [in] Roll-off of the filter
		const int Nfil=-1		//!< [in] Number of filter taps
		);


	//! Generate pilot signals
	void AddKnownSymbols();

	//! Returns the input block size of the transmitter in bits
	/*!
	Notice that this number depends on the active blocks*/
	int  GetInputSize() const
	{
		return Nbit;
	};

	//! Returns the out block size of the transmitter in samples. 
	/*! Notice that this number depends on the active blocks, in particular the presence of TX filter and of the pilots.
	*/
	int  GetOutputSize() const
	{
		return nsamp;
	};
    
	//! Returns the frequency spacing for CCCPM. 
	double GetDeltaf() const
	{
		return Deltaf;
	};

	//! Returns the code rate. 
	double GetRate() const
	{
		return rate;
	};

	
	//! Display the TX system parameters to the desired output.
	int Display(
		FILE* file=stdout //!< [out] Pointer to the output stream
		);


	//! Run the DVBRCS transmitter
	/** The sequence of TX processing blocks, as specified  during configuration, are applied to the input sequence of bits to produce the sequence of 
	complex points. The size of each input packet, as well as the size of the output packet depends on the active blocks.
	\return the number of generated samples \sa GetInputSize(), \sa GetOutputSize() */
	int Run(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		double *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		cmplx *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		)
	{
	return(Run(frames,data,(double *)tx));
	}	
#endif



	friend class DVBRCS_Receiver;

private:
	int CPM_flag;			//!< Specify if CCCPM or TCLM MODCOD is selected.
	int MODCOD;				//!< Current MODCOD

	TurboPHI_Encoder* Enc;  //!< Embedded TurboPHI encoder (only for TCLM)
	Modulator* LinMod;		//!< Embedded Linear Modulator (only for TCLM)
	Filter* TXFil;			//!< Embedded Shaping filter (only for TCLM)

	
	Interleaver*	Int;	//!< Embedded Interleaver (for both CCCPM and TCLM)

	Trellis_Encoder* TrEnc; //!< Embedded Convolutional Encoder (only for CCCPM)
	CPM_Modulator* CPMMod;	//!< Embedded CPM Modulator (only for CCCPM)
	Puncturer*   Punct;		//!< Embedded Puncturer (only for CCCPM)
	ShiftFrequency* shiftf; //!< Frequency shift for tilt removal (only for CCCPM)

	int period;				//!< Puncturing period
	int *pattern;			//!< Puncturing pattern

	int		*buffi1;		//!< Pointer to internal buffer storing  bits
	int		*buffi2;		//!< Pointer to internal buffer storing  bits
	double	*buffd1;		//!< Pointer to internal buffer storing constellation points
	double	*buffd2;		//!< Pointer to internal buffer storing constellation points

	int Nbit;				//!< Current number of input bits (before Encoder)
	int nsamp;				//!< Total number of transmitted samples
	int ns;					//!< Number of samples per symbol
	double rate;			//!< Code rate 
	int Nfil;				//!< Number of symbols in the filter (one side)
	double Deltaf;			//!< Carrier spacing normalized to the symbol rate

	int Nbps;				//!< Number of bits per symbol 
	int nsymb;				//!< Number of symbols
	
	int Nbit_int;			//!< Current number of coded bits (after Convolutional Encoder)
	int Nterm;				//!< Trellis termination symbols
	double alpha_RC;		//!< Phase response of the AV CPM pulse shape

	int rate_index;			//!< Index of Turbo-Phi code rate (0=1/3, 1=1/2, 2=2/3, 3=3/4, 4=4/5, 5=5/6, 6=6/7, 7=7/8)


	//! Set the active MODCOD of DVB-RCS transmitter for the TCLM interface
	int SetMODCOD_TCLM(
		const int MODCOD		//!< [in] Index of MODCOD
		);

	//! Set the active MODCOD of DVB-RCS transmitter for the CCCPM interface
	int SetMODCOD_CCCPM(
		const int MODCOD		//!< [in] ndex of MODCOD
		);
	
	//! Bit interleaving which defines the mapping between coded bits and modulation symbols
	void BitMapping(
		const int K,			//!!< [in] Number of information bits
		const int N,			//!!< [in] Number of coded bits
		const int m,			//!!< [in] Number of bits per symbol 
		const int rate_index,	//!!< [in] Index of code rate (0=1/3, 1=1/2, 2=2/3, 3=3/4, 4=4/5, 5=5/6, 6=6/7, 7=7/8)
		int* perm				//!!< [out] Output permutation
		);

	//! Retrieve the unique work for TCLM
	int GetUW(
		int MODCOD,				//!< [in] Index of MODCOD
		int *UW					//!< [out] Unique word
		);

	//! Run the DVBRCS transmitter for the CCCPM interface
	int Run_CCCPM(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		double *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		);

	//! Run the DVBRCS transmitter for the TCLM interface
	int Run_TCLM(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		double *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		);

#ifdef CTOPCOM
	//! Overload of Run_CCCPM with complex signals
		int Run_CCCPM(
			const int frames,	//!< [in] Number of transmitted frames
			const int *data,	//!< [in] Input data bits
			cmplx *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		)
		{
		return(Run_CCCPM(frames,data,(double*)tx));
		}

	//! Overload of Run_TCLM with complex signals
		int Run_TCLM(
			const int frames,	//!< [in] Number of transmitted frames
			const int *data,	//!< [in] Input data bits
			cmplx *tx			//!< [out] Transmitted complex signal (ns samples per symbol) 
		)
		{
		return(Run_TCLM(frames,data,(double*)tx));
		}
#endif

	bool withpil;			//!< Insert pilots?

	double *table_cos;
	double *table_sin;
	
	double rolloff;			//!< Roll-off factor of shaping filter
	
	int *UW;				//!< Unique word vector for TCLM
	int nUW;				//!< Size of UW
	
	int	Npre;				//!< Pre-amble length in symbols(for TCLM)
	int	Npost;				//!< Post-amble length ib symbols (for TCLM)
	int	Npil;				//!< Number of pilot symbols (for TCLM)
	int	Ppil;				//!< Period of pilot symbols

	
	static const int UW1=0x7CD593AD;	//!< First part of unique word for CCCPM
	static const int UW2=0xF7818AC8;	//!< Second part of unique word for CCCPM

	int *pilots;			//!< Binary vector containing pilot string
	
	// Static variables
	static int seedsc;
};
