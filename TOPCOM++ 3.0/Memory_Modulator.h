#pragma once
#include "Modulator.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class Memory_Modulator and its interfaces
*/
/*! \ingroup Modems
\brief Generic  modulator with memory

The class allows to define modulators with memory, where the transmitted constellation point at any given time
depends also on adjacent (left and right) \e L symbols.
This class is typically used in dynamic predistortion schemes to combat inter-symbol-interference.

The mehod SetParameters() initializes the block with a pointer to a reference Modulator (without memory)
and an integer specifying the memory \e L of the block.

The user can then modify the output by modifiying the elements of the public vector #constellation.

At each tic of the Run() method the \f$m\f$ input bits are converted into an 
M-ary integer and permuted through the law specified by the Modulator mapping. The resulting index is stored
in a length \f$ 2L+1 \f$ FIFO. The full content of the FIFO is then used to address the vector #constellation of size 
\f$ M^{2L+1} \f$ that provides the tranmitted constellation point.

\author Guido Montorsi
*/
class Memory_Modulator
{
public:
	Memory_Modulator(void);

	~Memory_Modulator(void);

	//! Run the modulator with memory for a given number of symbols.
	/*! If the input are in binary form (bin is true), for each symbol a block of \f$ m=\log_2 M \f$ consecutive bits are grouped to form 
	a label	of the symbol and then mapped to an index of modulation through the Modulator:mapping. 
	The resulting index is stored in a \f$ 2L+1 \f$ FIFO.  The full content of the FIFO is then used to address the buffer #constellation of size 
	\f$ M^{2L+1} \f$. When bin is false the input to the block is the sequence of labels 
	*/
	void Run(const int ntics,	//!< Number of generated symbols
		const int* inp,			//!< Pointer the input sequence (bits or integers in [0,M-1] depending on bin) 
		double* out,			//!< Pointer to sequence of constellation points
		const bool bin=true,	//!< Flag to indicate if input is binary
		int* ssymb=0,			//!< Optional output buffer to store sequence of supersymbol indexes
		int* symb=0			//!< Optional output buffer to store sequence of labels 
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,		//!< Number of processed input samples
		const int* inp,	    //!< Pointer the input sequence (bits or integers in [0,M-1] depending on bin)
		cmplx *out,			//!< Pointer to sequence of constellation points
		const bool bin=true,	//!< Flag to indicate if input is binary		
		int* ssymb=0,			//!< Optional output buffer to store sequence of supersymbol indexes
		int* symb=0			//!< Optional output buffer to store sequence of labels 
		)
	{
	Run(ntics,inp,(double*)out,bin,ssymb,symb);
	return;
	}
#endif

	//! Configure the block
	/*! The user passes a pointer to the memoryless Modulator and the parameter \e L that specifies
	the memory of the block.
	After initialization the block behaves like a memoryless modulator, with an additional latency of \e L symbols.
	The user can tune the \f$ M^{2L+1}\f$ constellation set by directly modifying the public buffer #constellation
	*/
	int SetParameters(
		const Modulator* mod,		//!< Pointer to the reference initial Modulator
		const int L					//!< Number of adjacent symbols included in memory
		);

	//! Pointer to the constellation of size \f$ M^{2L+1} \f$ used by the modulator
	/*! Modify the content of this buffer to update the constellation */
	double* constellation;


private:
	int* mapping;
	int* line;
	int time;
	int m;
	int L;
	int ML;
};
