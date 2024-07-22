// Viterbi.h: interface for the Viterbi class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Trellis.h"

/*! \file 
\brief Declaration of the class Viterbi
*/
/*! \ingroup Codecs
\brief Viterbi decoding for time invariant trellises. 

The parameters specified by the user are the pointer to a Trellis, 
the size of the total window or latency of
the Viterbi decoder \f$W\f$ and the updating window \f$U\f$, which have a default value of 10.

Each tic of the Run method  processes a single trellis step and thus requires at the input the metrics for the
\f$k\f$-step and output the decoded symbols at trellis step \f$k-W\f$.

In order to guarantee almost optimal results with the Viterbi decoder the difference 
\f$W-U\f$ must be greater than six to
seven time the constraint length of the trellis.

For an example of its use see e.g. the test program "test_Cont_Convolutional.cpp".

\author Guido Montorsi
*/
class Viterbi  
{
public:
	Viterbi(const Trellis* trellis,		//!< Pointer to a Trellis
		const int delay,				//!< Size of the total window or latency of the Viterbi decoder \f$W\f$
		const int ngroup=10				//!< Updating window size (backtracking period)
		);
	virtual ~Viterbi();

	//! Run the Viterbi decoder
	/** The user provide metrics on both input and output symbols to the encoder.

	The method return the ML input sequence and optionally (if the correspondent pointer is
	different from zero) the ML output sequence.
	 */
	void Run(const int tics,  //!< Number of trellis steps.
			const int *Inf,   //!< Input symbol metrics.
			const int *Cod,   //!< Coded symbol metrics
			int *OInf,	      //!< Decoded input symbols
			int* OCod=0,	  //!< Decoded output symbols
			bool block=false //!< Viterbi working in block mode.
			);

	//! Run the Viterbi decoder
	/** The user provide metrics only on output symbols to the encoder.

	The method return the ML input sequence and optionally (if the correspondent pointer is
	different from zero) the ML output sequence.
	 */
	void Run(const int tics,	//!< Number of trellis steps.
			const int* Cod,	//!< Coded symbol metrics
			int* OInf,		//!< Decoded input symbols
			int* OCod=0,		//!< Decoded Output symbols
			bool block=false	//!< Viterbi working in block mode.
			);
	friend class ML_Receiver;

private:
	int* forward;		// Starting state
	int* inputforw;		// Input symbol
	int* outputforw;	// Code Symbol


	int* buff_trace_back;	// trace back buffer
	int  win;				// Circular Window pointer

	int* dec;			// Decision buffer
	int  buffersize	;	// Decision buffer size

	int* napath;			// new path metrics
	int* apath;				// old path metrics

	int ns;				// Number of states
	int ni;				// Number of input symbols
	int no;				// Number of output symbols
	int ngroup;			// Decision window
	int window;			// Delay of Viterbi
};

