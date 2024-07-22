// RS_Decoder.h: interface for the RS_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include "RS_Encoder.h"
/*! \file 
\brief Declaration of class RS_Decoder.

*/
/*! \ingroup Codecs
\brief Reed- Solomon decoder.

	This class implements all functionalities of a RS hard decoder.
	The parameter of the RS code, which are derived from the corresponding 
	element of the class RS_Encoder, are the codeword length \f$n\f$,
	the information word length \f$K\f$ and the number of shortening symbols.

	For an example of its use see e.g. the test program "test_RS.cpp".


\author Gabriella Bosco
*/

class RS_Decoder  
{
public:
	RS_Decoder();
	virtual ~RS_Decoder();

	//! Set the main parameters of the code
	void SetParameters(
			RS_Encoder* encoder				//!< Member of the class RS_Encoder
			);
	//! Run the decoder
		void Run(
		const int blocks,	//!< Number of blocks of bits decoded. 
		const int* Input,	//!< Input symbols.
		int* Output		//!< Output symbols.
			);

private:
	RS_Encoder* code;	//!< Member of the class RS_Encoder containing the code parameters
	int nn,				//!< Codeword length.
	kk,					//!< Information word length.
	mm;					//!< Bits per symbol
    int ns;				//!< Number of shortening bits
	int* bufin;			//!< Buffer containing the vector of bits to be decoded.
	int *recd;			//!< Codeword
	int *lambda,		//!< Erasures Locator poly
		*s;				//!< Syndrome poly 
	int *b,				//!< Redundancy bits
		*t,				//!< Vector used in the Berlekamp-Massey algorithm
		*omega;			//!< Error+Erasures Locator poly
	int *root,			//!< Roots (index form)
		*reg,			//!< Vector used in the Chien search
		*loc;			//!< Error location numbers
};

