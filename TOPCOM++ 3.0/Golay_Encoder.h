// Golay_Encoder.h: interface for the Golay_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
/*! \file 
\brief Declaration of the class Golay_Encoder.

*/
/*! \ingroup Codecs
\brief Golay code (24,12) encoder.

For an example of its use see e.g. the test program "test_Golay.cpp".

\author Gabriella Bosco
*/

class Golay_Encoder  
{
public:
	Golay_Encoder();
	virtual ~Golay_Encoder();
	//! Run the encoder
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input signal. 
		int* Output			//!< Output signal.
		);
	int n,	//!< Codeword length.
		k;	//!< Information word length.
	friend class Golay_Decoder;
private:
long *encoding_table;  //!< Encoding table
};

