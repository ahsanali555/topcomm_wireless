// BCH_Decoder.h: interface for the BCH_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#include "BCH_Encoder.h"


#pragma once
/*! \file 
\brief Declaration of the class BCH_Decoder */
/*! \ingroup Codecs
\brief BCH decoder

The
parameters of the BCH code are derived from the corresponding
element of the class BCH_Encoder.  A single step of
the Run() method accepts a block of \f$n\f$ coded bits
 and produces a block of \f$k\f$ decoded bits. The flag withoutput 
 can be used to force the output of the block to be
the entire decoded codeword, composed by \f$n\f$ bits.

  	For an example of its use see e.g. the test program "test_BCH.cpp".

\author Gabriella Bosco
*/

class BCH_Decoder  
{
public:
	BCH_Decoder();  
	virtual ~BCH_Decoder(); 

	//! Set the main parameters of the code
	void SetParameters(
			BCH_Encoder* encoder	//!< Member of the class BCH_Encoder that descibes the code.
			);

	//! Run the decoder
	/*! The output buffer can coincide with the input buffer (if withoutput=true) */
		void Run(
			const int blocks,	//!< Number of encoded blocks 
			const int* Input,	//!< Input bits
			int* Output,		//!< Output bits
			const bool withoutput=false	//!< If true, the output is the decoded codeword instead of the decoded information vector. 
			);


		friend class PC_BCH_Decoder;
	BCH_Encoder* code;	//!< reference BCH Encoder

private:
	int n,				//!< Codeword length (before shortening).
	k,					//!< Information word length.
	t;					//!< Error correction capability of the code
	int nsh;			//!< Number of shortened bits
	int elp[1026][1026];//!< Error location polynomial (elp)
	int d[1026];		//!< Discrepancy at each step of the Berlekamp's algorithm
	int l[1026];		//!< Degree of the elp at each step of the Berlekamp's algorithm
	int u_lu[1026];		//!< Difference between	the step number and the degree of the elp
	int s[1026];		//!< Syndromes
	int root[200];		//!< Roots of the error location polynomial
	int loc[200];		//!< Error location number indices
	int reg[201];		//!< Vector used in the Chien search
	int *bufin,*bufout;	//!< Decoder input and output vectors 
};
