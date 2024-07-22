// Golay_Decoder.h: interface for the Golay_Decoder class.
//
//////////////////////////////////////////////////////////////////////
#include "Golay_Encoder.h"

#pragma once
/*! \file 
\brief Declaration of the class Golay_Decoder.
*/
/*! \ingroup Codecs
\brief Golay code (24,12) maximum likelihood decoder

 For an example of its use see e.g. the test program "test_Golay.cpp".

\author Gabriella Bosco
*/
class Golay_Decoder  
{
public:
	Golay_Decoder();
	virtual ~Golay_Decoder();

	//! Run the decoder
	void Run(
		const int	blocks,		//!< Number of encoded blocks
		const int*	Input,		//!< Input bits
		int* Output			//!< Output bits
		);
int n,	//!< Codeword length.
	k;	//!< Information word length.
private:
long *decoding_table;	//! Decoding table
int a[4];				//! Error pattern
};

