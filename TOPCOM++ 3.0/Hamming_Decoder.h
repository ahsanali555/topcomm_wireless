// Hamming_Decoder.h: interface for the Hamming_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#include "Hamming_Encoder.h"

#pragma once
/*! \file 
\brief Declaration of class Hamming_Encoder.
*/
/*! \ingroup Codecs
\brief Extended Hamming codes decoder.

For its use see the class PC_Decoder which has it as a member

\author Gabriella Bosco
*/

class Hamming_Decoder  
{
public:
	Hamming_Decoder();
	virtual ~Hamming_Decoder();

	//! Set the main parameters of the code
		void SetParameters(
			Hamming_Encoder* encoder,		//!< Element of the class Hamming_Encoder 
			bool siso=false				//!< If siso=true, the decoding table used by the SISO algorithm is generated
						);
	//! Run the decoder
		void Run(
		const int blocks,	//!< Number of blocks. 
		const int* Input,	//!< Input bits.
		int* Output		//!< Output bits.
			);
	
	//! SISO algorithm (to use it, the variable 'siso' in SetParameters must be 'true')	
		void RunSISO(
			const int blocks,	//!< Number of blocks. 
			const int* Input,	//!< Input sequence of LLRs
			int* Output		//!< Output sequence of LLRs
			);

		friend class PC_Decoder;
private:

		int n0,	//!< Codeword length.
		k0;	//!< Information word length.
		Hamming_Encoder* code;	//!< Member of the class Hamming_Encoder containing the code parameters
		int MAT[15000][4];	//!< Matrix containing the codewords of Hamming weight 4 having a 1 in the first position
		int nword;			//!< Number of codewords of Hamming weight 4 having a 1 in the first position
		int *bufin,*bufout,*buffer;
		int *s,*x,*r,*utmp;
		int *xtmp,*c;
		int *rel;
		int *DM;
		int *ytmp,*aff;
		int *indice;
		//!< Find the position of the minimum value in an integer vector v of length n
		int find_min(int *v,int n);
 
};


