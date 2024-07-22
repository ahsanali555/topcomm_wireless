// Reed_Muller_Encoder.h: interface for the Reed_Muller_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

/*! \file 
\brief Declaration of classes Reed_Muller_Encoder and Reed_Muller_Decoder.

*/
/*! \ingroup Codecs
\brief Generic Reed Muller Encoder

  Reed Muller encoders are defined through the two parameters  \e m and \e r that specify
  the base 2 logarithm of the code size and the "order" of the code.
  The free distance of the code is 2^(m-r).

For an example of its use see e.g. the test program "test_reed_muller.cpp".

\author Guido Montorsi
*/

class Reed_Muller_Encoder  
{
public:
	Reed_Muller_Encoder();

	//! Set the main parameters of code
	void SetParameters(
		const int m, //!< Log2 of code size
		const int r	 //!< Reed-Muller order
		);

	//! Run-time method
	void Run(const int tics,		//!< Number of encoded codewords
		const int* input,	//!< Input bits
		int*output			//!< Coded bits
		);
	virtual ~Reed_Muller_Encoder();
	int r;
	int m;	//!< Log2 of code length
	int n;	//!< Code length
	int k;	//!< Information length
	int *G;	//!< Generator matrix

};

/*! \ingroup Codecs
\brief Generic Reed Muller Decoder

Two versions of soft decoders are provided. The Run() method works only for order one encoders
and performs full decoding through Fast Hadamard transform. The RunSub() is a faster
suboptimal soft decoder that works decomposing 
the RM code as a generalized concatenated code (See the manual for more details).

For an example of its use see e.g. the test program "test_reed_muller.cpp".

\author Guido Montorsi

*/

class Reed_Muller_Decoder  
{
public:
	Reed_Muller_Decoder();

	//! Set the code parameters through a pointer to correspondent encoder.
	void SetParameters(
		const Reed_Muller_Encoder*  //!< Pointer to correspondent encoder.
		);

	//! Full Reed-Muller soft decoder (to be used only if r=1)
	void Run(const int tics,//!< Number of decoded codewords
		const int* llr,		//!< Input quantized LLR
		int 		*out	//!< Output decoded bits 
		);

	//! Reed-Muller soft decoder with suboptimal algorithm
	void RunSub(
		const int tics, //!< Number of decoded codewords
		const int* llr, //!< Input quantized LLR
		int *out		//!< Output decoded bits 
		);

	virtual ~Reed_Muller_Decoder();

private:
	const Reed_Muller_Encoder *cod;
	int* buff;//!< temporary buffer for Transform
	void GMC(const int r, const int m, const int *y, int*c);

};

