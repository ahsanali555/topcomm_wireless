// Block_Decoder.h: interface for the Block_Decoder class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "Block_Encoder.h"
/*! \file 
\brief Declaration of the class Block_Decoder.

*/
/*! \ingroup Codecs 
\brief Generic decoder for linear block codes

The Block_Decoder implements the decoding of a systematic block code using a 
 standard array decoding algorithm, based on the use of a decoding look-up table.
Since the decoding complexity grows exponentially with \f$n-k\f$, 
it is possible to simulate short codes only (\f$n<24\f$). Use ad-hoc algorithm if available.


 For an example of its use see e.g. the test program "test_Block_Codes.cpp".

  \author Gabriella Bosco
*/

class Block_Decoder  
{
public:
	Block_Decoder();
	virtual ~Block_Decoder();
	//! Set the main parameters of the code
	void SetParameters(const Block_Encoder* a,	//!< Encoder
							bool flagin=false  //!< If flagin=true, a correction is attempted only if the weight of the error is equal or lower than the correction capability 
				);
 void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input signal. 
		int* Output			//!< Output signal.
		);

	int n;		//!< Codeword length
	int k;		//!< Information vector length 
	int dmin;	//!< Minimum distance of the code

private:
	//! Evaluate the syndrome of the received vector 
	long get_syndrome(long pattern);   

	int t;			//!< Correction capability of the code
	long *matrix;   //!< Parity check matrix
	int *weigth;	//!< Haming weigth of error patterns
	bool flag;		//!< If flag=true, a correction is attempted only if the weight of the error is equal or lower than the correction capability 
	long *decoding_table;	//!< Decoding table (standard array)
};
