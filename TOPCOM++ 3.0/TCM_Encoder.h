/*! \file
\brief Declaration of TCM_Encoder class and its interface functions 

*/

#pragma once
#include "Trellis.h"
#include "Modulator.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file
\brief Declaration of class TCM_Encoder
*/ 

/*! \ingroup Modems 
\brief %Trellis Coded Modulation Encoder.

  The TCM encoder is specified through 3 basic constituent blocks:
  - A constellation of size \f$M \f$ specified through a pointer to the class Modulator
  - Each trellis section in the TCM encoder is associated to L
  modulation signals, hence the set of possible points is \f$M^L\f$
  - A set of kpar uncoded bits that generate parallel transitions
  - A trellis encoder that convolutional encode  k-kpar input bits into
  n output bits
  - An arbitrary mapper that assign a binary labeling to the set of modulation points

  This set of parameters allows to represents virtually all TCM encoder proposed in literature
  and have several degree of freedom. 
	
	For morre details see the manual.
  To get a generic "good" TCM modulator use the interface function GoodTCM(). 

  The encoder Run() method can work in continuous or terminated mode. 
  In continuous mode, each tic of the encoder
corresponds to a single step in the trellis diagram. At the end of each call the object preserves the state of the
encoder.

In terminated (block) mode, the encoder state is periodically (with period \f$N\f$) forced back to the
identity state. The forcing of the state requires to perform additional $n$ trellis steps, called
terminating steps, which are not associated to an input. As a result of the terminating operation
the rate of the encoder is slightly decreased and does not coincide anymore with the rate of the
original trellis diagram. The rate reduction obviously depends on the period $N$ and the number of
terminating steps $n$. In block mode a single tic of the encoder performs the encoding of a whole
block of data, and thus the encoder state is always zero when exiting from the Run() method.

If one is interested in block terminated encoding, the method SetTerminated() allows to specify
the size,  in terms of number of trellis steps, of the block. A flag indicates if the specified
size refers to the number of steps before the termination (\f$N\f$) or after the termination (\f$N+n\f$).
Terminated mode can not be used with the TCM\_Decoder.

For an example of its use see e.g. the test program "test_TCM.cpp" 

\author Guido Montorsi
*/
class TCM_Encoder  
{
public:
	TCM_Encoder();
	virtual ~TCM_Encoder();

	//! Initialize the TCM encoder 
	void SetParameters(
		Modulator* modin,		//!< Modulation used
		const int	nsig,		//!< Number of signals per trellis step
		Trellis * trel,			//!< Trellis decribing the encoder
		const int *map=0,		//!< Mapping of bits to modulation symbols
		const int kpar=0		//!< Number of uncoded bits (parallel transitions)	
		);
		
	//! Simulate the TCM encoder. 
	/*!
	  Warning: its behavior changes depending if the encoder is in terminated mode. 
	  (See TCM_Encoder::SetTerminated)
	  */
	 void Run(const int ntics,	//!< Number of tics
		 const int* inp,		//!< Input bits
		 double* out			//!< Outout samples
		 );

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,	//!< Number of tics
		 const int* inp,		//!< Input bits
		 cmplx* out				//!< Complex Outout samples
		 )
	{
	Run(ntics,inp,(double*) out);
	return;
	}
#endif


	 void Display(FILE* file=stdout);



	//! Force the encoder to work in terminated mode.
	/*!	Periodically the encoder is forced to go back to the zero state. 
	When operating in terminated mode the Run method change its behavior. The encoder
	in fact generates blocks of output signals. The number of input bits encoded
	is slightly decreased since terminating symbols are generated without any inputs.
	The rate of the encoder and its spectral efficiency are consequently slightly reduced.
	In terminated mode the TCM encoder work like a block encoder and the state is always zero
	at the beginning of encoding operations.
	*/
	void SetTerminated(
		const int nsymbols,	//!< Number of trellis steps for termination
		bool before=true	//!< Before or after termination ?
		);

	//! Reset the state of the encoder.
	void Reset();

	int K;						//!< input bits to the encoder
	int nsig;					//!< Number of signals per trellis step
	int nsteps;					//!< Number of trellis steps for termination
	bool terminated;
	friend TCM_Encoder* TCM_4D8PSK();
	int nsymbols;
	const int *map;				//!< Mapping
	friend class TCM_Decoder;
	friend class SCTCM_Encoder;
	friend class SCTCM_Decoder;
	friend class MIMO_Decoder;

private:
	Trellis* trel;				//!< Trellis describing the encoder
	Modulator* mod;				//!< Constellation
	int state;					//!< State of encoder
	int kpar;					//!< Uncoded bits (parallel transitions)
	int nterm;					//!< Number of terminating steps
	int* termination;



friend
TCM_Encoder* TCM_4D8PSK(int typ);
friend
TCM_Encoder* TCM_4D8PSK_5C(int type);
friend
TCM_Encoder* GoodTCM(const int k,const int L,const int mem, const int m,const int modtype);
friend
TCM_Encoder* Divsalar56();

};


/*! \ingroup Interface
\brief Return the TCM encoder specified by the standard CCSDS "ECSS-E-50-05A(24January2003)"

Differential precoding is not implemented. 
*/

TCM_Encoder* TCM_4D8PSK(
						int type		//!< (0): 2.0 b/s/Hz (1):2.5 bit/s/Hz
						);

/*! \ingroup Interface
\brief Return the TCM encoder specified by the standard CCSDS "ECSS-E-ST-50-05C_Rev.1(6March2009)"

Differential precoding is not implemented.
*/

TCM_Encoder* TCM_4D8PSK_5C(
						   int type		//!< (0): 2.0 b/s/Hz (1):2.5 bit/s/Hz
						);

/*! \ingroup Interface 
\brief Return  a "good" TCM encoder.

Table of good codes are taken from the literature.
*/

TCM_Encoder* GoodTCM(const int k,		//!< Input bits
					 const int L,		//!< Modulation signals.
					 const int mem,		//!< Memory of encoder
					 const int m,		//!< Number of bits of modulation signal.
					 const int modtype	//!< Modulation type
					 );
