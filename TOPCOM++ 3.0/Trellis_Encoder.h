// trellis_encoder.h: interface for the trellis_encoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Trellis.h"
#include <stdlib.h>

/*! \file 
\brief Declaration of the class Trellis_Encoder
*/
/*! \ingroup Codecs
\brief Generic encoder based on time-invariant trellis. 

This class implements an encoder whose input-output
relationships can be described through a \e time-invariant \e trellis diagram.

The initialization of the encoder is performed by passing the
pointer to the corresponding trellis. The
encoder Run() method works in continuous or terminated mode. In
continuous mode, each tic of the encoder corresponds to a single
step in the trellis diagram. At the end of each call the object
preserves the state of the encoder.

In terminated (block) mode, the encoder state is periodically (with period \f$N\f$) 
forced back to the identity state. The
forcing of the state requires to perform additional $n$ trellis steps, called terminating steps, 
which are not
associated to an input. As a result of the terminating operation the rate of the encoder is slightly decreased and does
not coincide  with the rate of the original trellis diagram. The rate reduction obviously depends on the period
\f$N\f$ and the number of terminating steps \f$n\f$. In block mode a single tic of the encoder performs the encoding of a whole
block of data, and thus the encoder state is always zero when exiting from the Run() method.


If one is interested in block terminated encoding, the method SetTerminated() allows to specify the size, in terms of
number of trellis steps, of the block. A flag indicates if the specified size refers to the number of steps before the
termination (\f$N\f$) or after the termination (\f$N+n\f$).

For an example of its use see e.g. the test program "test_Cont_Convolutional.cpp" 

\author Guido Montorsi
*/

class Trellis_Encoder  
{
public:
	Trellis_Encoder(Trellis*);
	virtual ~Trellis_Encoder();
	//! Run the encoder
	/*! The encoder Run method works in continuous or terminated mode. 
	
	  In continuous mode, each tic of the encoder corresponds to a single
step in the trellis diagram. At the end of each call the object
preserves the state of the encoder.

In terminated (block) mode, the encoder state is periodically 
(with period \f$N\f$) forced back to the identity state. The
forcing of the state requires to perform additional \f$n\f$ trellis steps, 
called terminating steps, which are not
associated to an input. 
As a result of the terminating operation the rate of the encoder is slightly decreased and does
not coincide  with the rate of the original trellis diagram. 
The rate reduction obviously depends on the period
\f$N\f$ and the number of terminating steps \f$n\f$. 
In block mode a single tic of the encoder performs the encoding of a whole
block of data, and thus the encoder state is always zero when exiting from the Run method.*/
	void Run(const int tics, //!< Number of tics
		const int* Input,	//!< Input symbols
		int* Output			//!< Output symbols
		);
	//! Specify the size of the block, in terms of number of trellis steps
	void SetTerminated(
		int tresteps,		//!< Number of trellis steps for termination 
		bool before=true	//!< Number of steps before termination or after termination ?
		);

	//! Set the encoder to work in tailbiting mode 
	/*! The first inf bits are used to determine the starting state that is forced to be
	identical to the final state. Works only for linear FF encoders */
	void SetTailbiting(
		int tresteps		//!< Number of trellis steps for tail-biting
		);

	int infbits;			//!< Number of information bits
	int codbits;			//!< Number of coded bits
	int nterm;
	Trellis* trel;
	int ttt[10];
private:
	int terminated;
	bool tailbiting;
	int* termination;
	int state;



/*! \ingroup Interface
\brief Return a binary Trellis_Encoder with a Trellis generated with the Canonical interface.
*/
friend
Trellis_Encoder *Binary_Convolutional_Encoder(int systematic, //!< Is the trellis systemtic?
				   const int k,		//!< Number of information bits
				   int n,			//!< Number of encoded bits 
				   const int* H,	//!< Feedback matrix
				   const int* Z		//!< Feedforward matrix
				   );
};


Trellis_Encoder *Binary_Convolutional_Encoder(int systematic, 
				   const int k,		
				   int n,			 
				   const int* H,	
				   const int* Z		
				   );