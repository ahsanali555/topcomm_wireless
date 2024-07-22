/*!\file
\brief Interface for the TCM_Decoder class.
*/

#pragma once
#include "TCM_Encoder.h"
#include "Delay.h"
#include "Demodulator.h"
#include "Viterbi.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of class TCM_Decoder
*/ 

/*! \ingroup Modems
\brief %Trellis Coded Modulation decoder.

  Implements a maximum likelihood decoder for trellis coded modulations
  specified through the class TCM_Encoder.
  Parameter of the block is the window of the embedded Viterbi decoder.
  The class also embeds the functionality of the demodulator and hence
  the input signals is the complex envelope of the received signal.

 
For an example of its use see e.g. the test program "test_TCM.cpp" 

\author Guido Montorsi
*/
class TCM_Decoder  
{
public:
	TCM_Decoder();
	virtual ~TCM_Decoder();

	//! Initialize the decoder.
	void SetParameters(const TCM_Encoder*, //!< Pointer to the correspondent TCM encoder 
		int window						   //!< Size of the window of the Viterbi Algorithm
		);

	//! Run the decoder
	void Run(const int ntics, //!< Number of trellis steps.
		const double* inp,	  //!< Sampled received signal (complex or real).
		int* out			  //!< Decoded bits.
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,	//!< Number of tics
		 const cmplx* inp,		//!< complex Sampled received signal
		 int* out				//!< Outout samples
		 )
	{
	Run(ntics,(double*)inp,out);
	return;
	}
#endif
private:
	Viterbi* vit;
	Demodulator* dem;
	const TCM_Encoder* enc;

	int *decpar;
	Delay *del;
	int *metrics;
	int* lf;
	int* signl;
};
