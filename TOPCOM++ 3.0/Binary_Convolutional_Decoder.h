// Binary_Convolutional_Decoder.h: interface for the Binary_Convolutional_Decoder class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include "Viterbi.h"
#include "Trellis.h"

/*! \file 
\brief Declaration of the class Binary_Convolutional_Decoder.*/
/*! \ingroup Codecs
\brief Convolutional decoder for \e binary codes.


The parameters specified by the user are the pointer to a Trellis, the size of the total window or latency of
the Viterbi decoder \f$W\f$ and the updating window \f$U\f$, which have a default value of 10.

Each tic of the Run method  processes a single trellis step and thus requires at the input the metrics for the
\f$k\f$-step and output the decoded bits at trellis step \f$k-W\f$.

 For an example of its use see e.g. the test program "test_Convolutional.cpp".

\author Guido Montorsi
*/

class Binary_Convolutional_Decoder  
{
public:
	Binary_Convolutional_Decoder();
	virtual ~Binary_Convolutional_Decoder();
	
	//! Set the parameters of the decoder
	void SetParameters(
		const Trellis* trellis,	//!< Pointer to a Trellis
		const int delay,		//!< Total latency of Viterbi decoder
		const int ngroup=10,	//!< Updating window size (backtracking period)
		bool blockin=false,		//!< Block or continuous  decoding ? 
		int nblock=0,			//!< Total number of information symbols in a block 
		int termination=0		//!< Number of symbols used for block termination
						);
    
	//! Run the decoder
	void Run(
			const int ntics,			//!< Number of trellis steps
			const int* Input,			//!< Input bits
			int* Output					//!< Decoded bits
				  );
    private:
	Viterbi	*Decoder;
	int *split;
	bool block;
	int nb_out,nb_in;
	int n;
	int nu;				//<! Number of termination symbols
	int *bufin;
	int *bufout;
};
