// PC_Encoder.h: interface for the PC_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#include "Hamming_Encoder.h"

#pragma once

/*! \file 
\brief Declaration of class PC_Encoder
*/
/*! \ingroup Codecs
\brief Generic product-code based on single parity and extended Hamming codes (enhanced TPC as option).

Product Codes are composed by the concatenation of two extended Hamming codes or parity check codes 
\f$C^1\f$ and \f$C^2\f$ with
parameters \f$(n_1,k_1)\f$ and \f$(n_2,k_2)\f$ respectively. The product code \f$P\f$=\f$C^1\times C^2\f$ is obtained by:
-# placing \f$(k_1\times k_2)\f$ information bits in an array of \f$k_1\f$ rows and \f$k_2\f$ columns,
-# coding the \f$k_1\f$ rows with code \f$C^2\f$,
-# coding the \f$k_2\f$ columns with code \f$C^1\f$.
The encoding of single constituent codes of type extended Hamming code is performed by the 
Hamming_Encoder class.
-# [Enhanced] A third  set of parity checks is optionally added by taking the modulo-2 sum of all 
diagonals. In this case the column code is shortened by one, so as to keep the total length to \f$N=n_1 n_2\f$.
The input block size is reduced to \f$k_1 (k_2-1)\f$.


The user specifies:
- The codeword lengths  \f$n_1\f$ of the first constituent code, which has to be
a power of 2 in the range (8,256).
- The codeword lengths  \f$n_2\f$ of the second constituent code, which has to be
a power of 2 in the range (8,256).

- Two additional optional flags are used to indicate if the constituent
code is a parity check code or an extended Hamming code (default)
- A final flag specifies if the TPC is "enhanced". In enhanced TPC,
the column code is shortened by one and the final row is composed by the parity checks of all diagonals of the matrix. 
The code size then remains the same but the rate is decreased.

For an example of its use see e.g. the test program "test_PC.cpp".

\author Gabriella Bosco
\author Guido Montorsi

*/

class PC_Encoder  
{
public:
	PC_Encoder();
	virtual ~PC_Encoder();

	//! Set the main parameters of the code
	void SetParameters(
			const int n1,	//!< First codeword length (must be equal to \f$2^m\f$)
			const int n2,    //!< Second codeword length (must be equal to \f$2^m\f$)
			const bool p1=false, 	//!< Flag to set a parity check code for the row code (def:extended Hamming code)
			const bool p2=false, 	//!< Flag to set a parity check code for the column code (def:extended Hamming code)
			const bool enh=false 	//!< Flag to set enhanced mode for TPC)
			);
	
	//! Run the encoder
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input signal. 
		int* Output			//!< Output signal.
		);

	friend class PC_Decoder;
	
	int n1,n2,	//!< Codeword length.
		k1,k2;	//!< Information word length.
	bool enh;
private:
	bool p1,p2;
	Hamming_Encoder *enc1,*enc2;	//!< Members of class Hamming_Encoder defining the two constituent code
	int *bufin,						//!< Input vector to the Hamming encoder
		*bufout;					//!< Output vector from the Hamming encoder
	int *prod;						//!< Product matrix
	int *vtmp;						//!< Internal buffer
};

//! Block Turbo Code (BTC) of Wimax 802.16e standard
PC_Encoder* BTC_WIMAX(
					  const int rate		//!< rate (12, 34, 35, 45, 23, 56)
					  );


