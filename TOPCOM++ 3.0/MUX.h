// MUX.h: interface for the MUX class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of the class MUX
*/
/*! \ingroup DSP 
\brief Multiplexer or generic Parallel to Parallel converter


The input samples are stored in an internal  buffer until \f$M\f$ of
them are collected. The content of the buffer is then written to the
output.

The block is initialized by the constructor with the parameter \f$M\f$
that denotes the size of output blocks.

The Run() method requires to specify the number of input samples,
and the pointers to the input and output buffer. It returns the
number of written output blocks.

The SetOffset() method allows to specify the initial offset of input stream,
which by default is zero.

\author Guido Montorsi
*/
class MUX  
{
public:
	//! Construction and setting of the output size
	MUX(const int sizeout //!< Size of the output block
		);

	//! Run the MUX on integers. Returns the number of output blocks
	int Run(
		const int nsamp,		//!< Number of input samples 
		const int*inp,			//!< input samples
		int*out					//!< output blocks
		);

	//! Run the MUX on doubles. Returns the number of output blocks
	int Run(
		const int nsamp,			//!< Number of input samples 
		const double *inp,			//!< input samples
		double *out					//!< output blocks
		);

	//! Set the initial offset of output blocks
	void SetOffset(const int off){time=off%sizeout;}
	virtual ~MUX();

#ifdef CTOPCOM
	int Run(
		const int nsamp,			//!< Number of input samples 
		const cmplx *inp,			//!< input samples
		cmplx *out					//!< output blocks
		)
	{
		return Run(
			nsamp,			//!< Number of input samples 
			(const double *)inp,			//!< input samples
			(double *)out					//!< output blocks
			);
	}
#endif
private:
	int time;
	void* buffer;
	int sizeout;
};

