#pragma once
#include "LDPC_Encoder.h"
/**************************** Derived ************************/
/*! \file 
\brief Declaration of classes LDPC_Encoder_CCSDS
*/
/*! \ingroup Codecs*/
/*@{*/
/*!
\brief Low Density Parity Check Encoder for CCSDS Deep space codes 131.1-o-1 

	The class LDPC_Encoder_CCSDS, derived from the class
	LDPC_Encoder, has been specifically developed for implementing
	the LDPC encoder described in this standard. The encoders of this
	document in fact do not fall in the generic description of a
	LDPC_Encoder, which requires the parity check matrix to be in
	lower triangular form.

	The SetParameters() method of this class only requires to specify
	the rate and the length of the code. 

	The Run() method is different from that in the base class, which
	proceeds sequentially by back substitution. The encoder algorithm is
	described in the standard, section 3.4.

	The decoder is the generic LDPC Decoder.

	For an example of its use see e.g. the test program "test_LDPC.cpp".

\author Guido Montorsi
*/

class LDPC_Encoder_CCSDS:public LDPC_Encoder
{

public:
	
	~LDPC_Encoder_CCSDS();
	void SetParameters(const int rate,		//!< Rate of code (12,23,45)
						const int length	//!< Information block length (1024, 4096, 16384)
					  );
	void Run(const int ntics,	//!< Number of encoded blocks
				const int* inp, //!< Input bits
				int* out		//! Output bits
				);
private:
	int *p1,*p2;
	int M;		//!< Size of subblocks
	int rate;	//!< Rate of code
	static const int par[];
	static const int par2[];
	static const int H[];
	static const int H2[];
	int *P[26];	// Permutations

};

