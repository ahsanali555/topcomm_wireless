// FIFO.h: interface for the FIFO class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class FIFO
*/
/*! \ingroup DSP 
\brief First Input First Output line (asynchronous buffer)

A FIFO buffer is required
between two blocks that provide or accept data at variable
instantaneous rate but with equal average rate. In formulae, the
two sequences of data block sizes \f$O_j\f$ and \f$I_j\f$  generated and
accepted by the two connected blocks must satisfy
\f[
\sum_{j=-\infty}^n I_j =\sum_{j=-\infty}^n O_j + A_n \;\;\forall
n\;\; |A_n|\leq A<\infty
\f]

In this case, a FIFO with size \f$2A\f$ is sufficient to guarantee
correct data flowing between the two blocks.

The FIFO is realized with a circular buffer of size \f$N\f$ where the
pointer of the input data and output data ($I$ and $O$) are
initially displaced by an amount \f$a\f$, e.g. \f$I=0\f$ and \f$O=a\f$.  Each
call to the Run() method specifies the input data size \f$NI\f$ and
the output data size $NO$. The input data and output data are
sequentially (for all \f$i=0,\ldots,\max(NI,NO)-1\f$) stored into the
buffer at position \f$(I+i)\bmod N\f$ and read out from the buffer at
positions \f$(O+i)\bmod N\f$.

If \f$i\geq NI\f$ the data are only read out. If \f$i\geq NO\f$ the data are
only stored.

At the end the pointer are updated as
\f[ I := (I+ NI)\bmod N\f]
\f[ O := (O+ NO)\bmod N \f]

  \author Guido Montorsi
*/
class FIFO  
{
public:
	FIFO();
	//! Set the main parameters
	void SetParameters(
		const int size,		//!< Size of buffer
		const int delay=-1	//!< delay of output stream (default to size/2)
		); 

	//! Run the FIFO with a prescribed number of input samples and output samples
	int Run(const int sampi,		//!< Number of provided input samples or blocks
		const int* input,			//!< Input samples
		const int sampo,			//!< Number of required output samples or blocks
		int* output					//!< Output samples
		);

	//! Run the FIFO with a prescribed number of input samples and output samples
	int Run(const int sampi,		//!< Number of provided input samples or blocks
		const double* input,			//!< Input samples
		const int sampo,			//!< Number of required output samples or blocks
		double* output					//!< Output samples
		);

	//! Return the current FIFO depth
	int Depth();

	virtual ~FIFO();
	int pui;
	int puo;

#ifdef CTOPCOM
	int Run(const int sampi,		//!< Number of provided input samples or blocks
		const cmplx* input,			//!< Input samples
		const int sampo,			//!< Number of required output samples or blocks
		cmplx* output					//!< Output samples
		)
	{
		return Run(sampi,		//!< Number of provided input samples or blocks
			(const double*) input,			//!< Input samples
			sampo,			//!< Number of required output samples or blocks
			(double*) output					//!< Output samples
			);
	}
#endif

private:
	int NI;
	int NO;
	void *buff;
	int size;
	int flag;
};

