// PC_Decoder.h: interface for the PC_BCH_Decoder class.
//
//////////////////////////////////////////////////////////////////////
#include "PC_Encoder.h"
#include "Hamming_Decoder.h"
#pragma once

/*! \file 
\brief Declaration of class PC_Decoder
*/
/*! \ingroup Codecs
\author Gabriella Bosco
\brief Iterative soft decoder for product-codes based on extended Hamming codes.


The user specifies:
- The number of iteration of the decoding algorithm. 
- The values of the feedback coefficients \f$\alpha[n]\f$. If they are not specified, 
all the coefficients are set to 0.5.

\see PC_Encoder for a description of the encoding scheme

	For an example of its use see e.g. the test program "test_PC.cpp".

  \author Guido Montorsi
*/


class PC_Decoder  
{
public:
	PC_Decoder();
	virtual ~PC_Decoder();

	//! Set the main parameters of the decoder
	void SetParameters(
			PC_Encoder* encoder,	//!< Pointer to the reference PC_Encoder
			int nit,				//!< Number of iterations in the decoding algorithm
			double *gains=0			//!< Vector containing the gain coefficients used in the decoding algorithm (default: all 0.5)
			);

	//! Run the decoder
	/*! The optional final pointer allows to specify the transmitted codeword. When present, the 
	iterative decoder is stopped if no errors are found. Otherwise all iterations are always performed.
	*/
	int Run(
		const int blocks,	//!< Number of decoded blocks
		const int* Input,	//!< Input signal
		int* Output,		//!< Output signal
		const int *data = 0	//!< Encoded bits for genie-aided stopping rule [optional]
			);
			


private:

		Hamming_Decoder *dec1,*dec2;	
		int n1,n2,			//!< Codeword length of constituent codes.
		k1,k2;				//!< Information word length of constituent codes..
		PC_Encoder* code;	//!< Member of the class PC_Encoder containing the code parameters
		int Niter;			//!< Number of iterations in the decoding algorithm
		double *AA;			//!< Vector containing the gain coefficients used in the decoding algorithm
		int *bufin,			//!< Input vector to the Hamming decoder
			*bufout;		//!< Output vector from the Hamming decoder
		int *R,			//!< llrs at eah iteration
			*W1,			//!< Extrinsic information matrix
			*W2,			//!< Extrinsic information matrix
			*W3;			//!< Extrinsic information matrix

		int Soft_Parity_Check(
			const int n, 
			const int *inp, 
			int *out, 
			const int beta=3);
		
};
