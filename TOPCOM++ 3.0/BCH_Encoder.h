// BCH_Encoder.h: interface for the BCH_Encoder class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <stdio.h>
/*! \file 
\brief Declaration of the class BCH_Encoder*/
/*! \ingroup Codecs 
\brief  BCH encoder

The BCH_Encoder class implements a programmable encoder for Bose-Chaudhuri-Hocquenghem (BCH)
codes. 
The user specifies:
- The codeword length n, which must be equal to \f$2^m-1\f$, where \f$m\leq 16\f$ is a positive integer.
- The correction capability of the code \f$t\f$ in number of bits.

  	For an example of its use see e.g. the test program "test_BCH.cpp".

\author Gabriella Bosco
*/


class BCH_Encoder  
{
public:
	BCH_Encoder(); 
	virtual ~BCH_Encoder(); 
	
	//! Set the main parameters of the code
	void SetParameters(
			int n,				//!< Codeword length (before shortening, must be equal to \f$2^m-1\f$)
			int t,				//!< Correction capabilty of the code
			int nsh=0			//!< Shortening bits
			);

	//! Run the encoder
	/*! The parity check bits are inserted before the systematic bits */
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input bits. 
		int* Output				//!< Output bits.
								/*!< The parity check bits are inserted before the systematic bits */
		);
	
	//! Display generator polynomial 
	void Display_Gen(FILE * file=stdout) const;

	friend class BCH_Decoder;
	
	int n,	//!< Codeword length.
		k,	//!< Information word length.
		t,	//!< Error correction capability of the code
		nsh;//!< Shortening bits
private:
	void read_p(int m);						//!< Read the generator polynomial of the Galois field \f$GF(2^m)\f$
    void generate_gf(int m);				//!<Generate field \f$GF(2^m)\f$ from the irreducible polynomial
    void gen_poly(int m,int length,int t);  //!<Compute the generator polynomial of a binary BCH code
	int *p;					//!< Galois field generator polynomial
	int *alpha_to;			//!< Look-up table for the conversione of elements of \f$GF(2^m)\f$ from index to polynomial form
	int *index_of;			//!< Look-up table for the conversione of elements of \f$GF(2^m)\f$ from polynomial to index form
	int *g;					//!< Code generator polynomial
	int *bufin,*bufout;		//!< Encoder input and output vectors 
};
