// PCCC_Decoder.h: interface for the PCCC_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include "PCCC_Encoder.h"
#include "SISO_Decoder.h"
/*! \file 
\brief Declaration of class PCCC_Encoder
*/
/*! \ingroup Codecs
\brief Generic binary Parallel Concatenated Convolutional Code Encoder.

	The PCCC_Decoder class is associated to the correspondent
	PCCC_Encoder.

	Through the method SetParameters the user specifies:
	- A pointer to the correspondent PCCC_Encoder that describes the code
	- The number of iterations \f$n_I\f$ of the decoder
	- The size of updating and training window \f$T_W\f$ and \f$U_W\f$ of the embedded SISO decoders. 
	(Optional, default to automatic computation)
	- The scaling factor \f$f\f$ of the input quantized LLRs. (Optional, default to 8, 
	which gives almost ideal performances).

	Optionally, with the method SetStop(), the user can activate or deactivate a stopping rule of the decoder. 

	The Run method takes as input the quantized LLRs of the coded bits 
	after a possible depuncturing, which inserts zero values for those bits that were not transmitted, 
	the iterative decoder
	is performed for \f$n_I\f$ iterations. 
	The output buffer is then filled with the last computed  LLRs of the information
	bits.

	For an example of its use see e.g. the test program "test_PCCC.cpp".

\author Guido Montorsi
*/
class PCCC_Decoder  
{
public:
	PCCC_Decoder();
	virtual ~PCCC_Decoder();

	//! Set the main parameters of the decoder
	void SetParameters(const PCCC_Encoder* refcodin,	//!< Reference PCCC encoder
						const int niter,			//!< Number of iterations
						const int nini=-1,				//!< Training window for SISOs (-1=automatic computation)
						const int ngroup=-1,			//!< Updating window for SISOs (-1=automatic computation)
						const double factor=8.			//!< factor for max* operation
						);

	//! Activate or deactivate a stopping rule for the decoder
	/*! The stopping rule is based on the comparison of the LLRs signs at the input and at the output of the lower
		SISO processor. If the signs agree for all bits the iterative decoding is stopped.
		Alternatively, if the user specifies the transmitted data to the Run method, the stopping
		rule is genie aided.*/
	void SetStop(const bool stopin=true){stop=stopin;}

	/*! Change the number of itreations of the decoder */
	void TuneParameters(const int niter);

	//! Run the decoder, one tics decode a whole codeword
	int Run(const int tics,			//!< Number of decoded codewords
		const int* inp,					//!< Input LLR of coded bits
		int* out,						//!< Output LLR of information bits
		const int* data=0,				//!< Optional transmitted data for genie-aided stopping rule
		const bool reset=true				//!< Optional transmitted data for genie-aided stopping rule
		);

	//! Run the depuncturer of the decoder 
	/*! With Complementary Code Combining, multiple  codewords associated to the
	same information blocks are transmitted through different channels with possibly different rates.
	The LLR of the bits are all summed up. The decoder run only when all LLR have been combined (added).
	This method runs only the depuncturer in front of the PCCC decoder and return the buffer with the 
	combined LLR
	*/
	void RunDepunct(const int tics,		//!< Number of depunctured codeword
		const int* llr,					//!< Input LLR of punctured coded bits
		int* out						//!< Output LLR of accumulated depunctured coded bits
		);
	//! Turn on and off the depuncturing algorithm
	/*! If the encoder has a puncturer associated with it, the depuncturer 
	algorithm is automatically turned on.
	*/
	bool puncturing;
	int* apriori; //!< Pointer to a-priori LLR on inf bits 

private:
	friend class RX_BGAN_102744;
	bool CCC;
	const PCCC_Encoder *refcod;		//!< Reference to the correspondent encoder.
	int niter;
	SISO_Decoder * uSISO;
	SISO_Decoder * lSISO;
	int *bufllr;		// Buffer to store depunctured llr;
	int *ext1;
	int *ext2;
	int ncod;
	bool stop;
	int sizebuf;
};

