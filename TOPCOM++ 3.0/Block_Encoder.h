// Block_Encoder.h: interface for the Block_Encoder class.
//
//////////////////////////////////////////////////////////////////////


#pragma once

/*! \file
\brief Declaration of the class Block_Decoder.
 */

 /*! \ingroup Codecs 
\brief  Generic linear block-codes encoder (systematic codes)

The Block_Encoder class implements the encoding of a systematic block code 
defined through its parity-check matrix 
\f$
{\bf H}=\left[{\bf P}|{\bf I}\right]
\f$
where \f${\bf I}\f$ is the \f$(n-k)\times(n-k)\f$ identity matrix and 
\f${\bf P}\f$ is a \f$(n-k)\times(k)\f$ binary matrix.

The user specifies:
- The codeword length \f$n<24\f$.
- The information vector length \f$k\f$
- The name of a file containing the non diagonal part \f${\bf P}\f$ 
of the parity check matrix of the code ((\f$n-k\f$) rows and \f$k\f$ columns).

The encoding is based on the use of a look-up table.

For an example of its use see e.g. the test program "test_Block_Codes.cpp".


\author Gabriella Bosco

*/
class Block_Encoder  
{
public:
	Block_Encoder();
	virtual ~Block_Encoder();
	
	//! Set the main parameters of the code
	void SetParameters(int nin,				//!< Codeword length (no greater than 24)
						int kin,			//!< Information vector length 
						char * namefile	//!< Name of the file containing the non diagonal part of the parity check matrix of the code ((n-k) rows and k columns)
						);
	//! Run the encoder
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input bits. 
		int* Output				//!< Output bits.
		);

	int n;		//!< Codeword length 
	int k;		//!< Information vector length 
	int dmin;	//!< Minimum distance of the code

	friend class Block_Decoder;

private:
	long *encoding_table;	//!< Look-up table for encoding
	long *matrix;			//!< Parity-check matrix

};


