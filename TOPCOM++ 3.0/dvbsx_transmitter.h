#pragma once
#include "PN_Source.h"
#include "BCH_Encoder.h"
#include "LDPC_Encoder.h"
#include "Interleaver.h"
#include "Modulator.h"
#include "Filter.h"
#include "Transponder.h"
#include "Reed_Muller_Encoder.h"

#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of the class DVBSX_Transmitter
*/

/*! \ingroup Systems
\brief Transmitter system as specified in the standard ETSI EN 302 307 V1.3.1 (2013-03)
with the extensions introduced in DVB-SX. 

The user is assumed to be familiar with the standard.
The TX system includes the following blocks:

- [optional] Energy dispersal 
- [optional] BCH Encoding 
- [optional] LDPC Encoding and Interleaving 

-  Mapping of bits to DVB-S2/SX constellation sets 

- [optional] Physical Layer Header Insertion
- [optional] Pilots Insertion
- [optional] PL scrambling
- [optional] Shaping with Square Root Raised Cosine filter of variable roll-off


All optional blocks can be inserted or removed toggling the correspondent flags 
or using the methods for their configuration. 
Notice  however that if BCH coding is inserted then also LDPC is inserted independently from the correspondent flag.

The Physiscal layer scrambler is configured with the method AddPLScrambler().
The square root raised cosine transmitting filter is configured with the method AddTXFilter().

The configuration of the transmitter is set with the method SetMODCOD(). 
All the previously described configuration flags must be called before the call to this method.

The transmitter also support static data predistortion (TunePredistortion())

\sa DVBSX_Receiver

\author Guido Montorsi
*/
class DVBSX_Transmitter
{
public:
	DVBSX_Transmitter(void);
	~DVBSX_Transmitter(void);

	bool withbitscrambler;		//!< Insert bit scrambler? (energy dispersal)
	bool withBCH;				//!< Insert BCH?
	bool withLDPC;				//!< Insert LDPC?
	bool withPLheader;			//!< Insert Physical Layer Header?
	bool withpil;				//!< Insert pilots?


	//! Add and configure the physical layer scrambler.
	/*! 
	Set n=0 to remove the physical layer scrambler.	*/
	void AddPLScrambler(
		const int n=1	//!< [in] Seed of the Gold sequence for scrambling
		);

	//! Add and configure the Square Root Raised Cosine shaping filter
	/*! 
		Set ns =1 to remove the filter. */
	int AddTXFilter(const int ns,	//!< [in] Number of samples per symbol (ns=1 means no filter)
		const double rolloff,		//!< [in] Roll-off of the SRRC filter
		const int Nfil=-1
		);

	//! Set the active MODCOD of DVB-SX transmitter
	//! \return an integer that is zero if the MODCOD has not been set.
	//! It must be called after the configuration of the active blocks in the transmitter.
	int SetMODCOD(
		const int MODCOD,		//!< [in] MODCOD Index (S2:1-28, SX:129-248) 
								//!< 129 and 131 specify very low SNR MODCODs sets and require the following index
		const int VLSNRMODCOD=0 //!< [in] Optional index to choose VLSNR MODCODS (129: 0-5 or 131: 0-3) 
		);



	//! Returns the reference Es/N0 in [dB] for Quasi Error Free over Linear channel. The reference values are those listed in the standard,.
	double GetEsN0_Lmin() const
	{if(index>=0)return refesn0L[index];else return 0.;} 

	//! Returns the reference Es/N0 in [dB] for  Quasi Error Free over Hard Limiter Non linear channel. 
	/*! available only for SX MODCODs (MODCOD>128). */
	double GetEsN0_NLmin() const
	{if(index>=0)return refesn0NL[index];else return 0.;}

	
	//! Display the TX system parameters to the desired output.
	void Display(FILE* file=stdout //!< [out] Pointer to the output stream
		) const;

	//! Tune the static modulation predistortion.
	/*!
	The predistortion modifies the position of transmitted constellation points so that 
	the centroids at the output of the non linear channel, which is
	supplied through a pointer to a Transponder, are a scaled version of 
	the nominal constellation points. The scaling value is specified by the parameter OBO.
	*/
	void  TunePredistortion(
		const int nsamples, //!< [in] Number of samples for tuning the predistortion
		Transponder* TTT,	//!< [in] Pointer to the used transponder (IMUX+TWT+OMUX)
		const double OBO,	//!< [in] Target scaling of output constellation [dB]
		const int off=-1	//!< [in] Delay of satellite chain in samples (-1: unknown)
							/**< If the delay is left unspecified (-1) the block evaluates internally
							the delay  by computing the peak of the total discrete time impulse response*/
		);

	//! Return the total back to back delay in samples, including nonlinearity.
	/*! The method evaluates and return the back to back delay of a transmission system embedding the DVB-SX transmitter, a Transponder supplied by the user,
	and the receiver. The delay is obtained  computing the peak of the total discrete time impulse response. 
	\return the back to back delay of a transmission system embedding the DVB-SX transmitter.	
	*/
	int  GetDelay(Transponder* TTT  //!< [in] Pointer to the transponder included in the system for evaluating the back to back delay.
									//!<
		) const;

	//! Run the DVB-SX transmitter
	/** The sequence of TX processing blocks, as specified  during configuration, are applied to the input sequence of bits to produce the sequence of complex points.
	The size of each input packet, as well as the size of the output packet depends on the active blocks.
	\return the number of generated samples 
	\sa GetInputSize(), \sa GetOutputSize() */
	int Run(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		double *tx,			//!< [out] Transmitted complex signal (ns samples per symbol 
		int *ref=0			//!< [in] Indexes of transmitted constellation points (used for tuning the receiver)
		);

	int* dataldpc;			//!< buffer to store output data of ldpc (used to stop iteration at the decoder)
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(
		const int frames,	//!< [in] Number of transmitted frames
		const int *data,	//!< [in] Input data bits
		cmplx *tx,			//!< [out] Transmitted complex signal (ns samples per symbol 
		int *ref=0			//!< [in] Index of transmitted constellation
		)
	{
	return(Run(frames,data,(double*)tx,ref));
	}	
#endif

	//! Returns the input block size of the transmitter in bits
	/*!
	Notice that this number depends on the active blocks*/
	int  GetInputSize() const
	{
		if(withBCH)			return K;
		else if(withLDPC)	return K1;
		else				return Ns*Mod->m;
	};

	//! Returns the out block size of the transmitter in samples. 
	/*! Notice that this number depends on the active blocks, in particular the presence of TX filter and of the pilots.
	*/
	int  GetOutputSize() const
	{
		return nsamp;
	};

	int Ns;					//!< Number of constellation data symbols 
	Modulator* Mod;			//!< Embedded Modulator
	Modulator* predMod;		//!< Embedded predistored Modulator
	int VLSNR;		//!< Index VLSNR

private:
	int MODCOD;					//!< Current MODCOD
	int VLSNRMODCOD;			//!< Current Very-Low SNR MODCOD

	//! Return the length and the values of the sequence of pilots
	/*! A value zero means no pilot in that position. */
	int GetPilot(
		int &nsymb,			//!< [out] Reference to the integer where the length of pilot sequence will be stored 
		double* pilot		//!< [out] Pointer to the buffer where pilot sequence(complex) will be stored 
		) const;


							/**< It is used \f$ M \f$ */
	bool SF;		//!< Flag to indicate spreading factor 2 on modulation (Bits are repeated before modulation)
	int K;			//!< Number of input bits (before BCH)
	int K1;			//!< Number of input bits (before LDPC)
	int N;			//!< Number of coded bits (after LDPC)
	int nsymb;		//!< Number of symbols (including pilots and header)
	int nsamp;		//!< Number of outputs samples 
	int S2;			//!< Specify if S2 or SX MODCOD is selected.

	BCH_Encoder*	BCH;	//!< Embedded BCH encoder
	LDPC_Encoder*	Enc;	//!< Embedded LDPC encoder	
	Interleaver*	Int;	//!< Embedded Interleaver
	Filter* TXFil;			//!< Embedded Shaping filter
	Reed_Muller_Encoder* PLS;	//!< Embedded Encoder 

	friend class DVBSX_Receiver;
	int index;				//!< Sequential index of MODCOD	
//	bool Normal;			//!< Normal or short SX mode	
	int		*buffi1;		//!< Pointer to internal buffer storing  bits
	int		*buffi2;		//!< Pointer to internal buffer storing  bits
	double	*buffd1;			//!< Pointer to internal buffer storing constellation points
	double	*buffd2;			//!< Pointer to internal buffer storing constellation points
	//! Physical layer scrambler
	int xi;			//! < s
	int yi;
	void PLscramble(int tics, const double* inp, double* out, bool inv=false) const;
	int nsc;	//!< Number of Gold sequence for PL scambler
	int ns;
	int rate;
	int blen;
	double rolloff;

	// Static variables
	static const int seedsc;
	static const int colord[];
	static const double refesn0L[];
	static const double refesn0NL[];
	static const int SOF;
};
