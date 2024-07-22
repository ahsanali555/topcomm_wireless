// BER_meter1.h: interface for the BER_meter class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include "Delay.h"
#include "Include.h"
#include <time.h>
#include <stdio.h>

/*! \file 
\brief Declaration of the class BER_meter
*/
/*! \ingroup Measurements 
\brief  Measure Bit Error Rate and Frame Error Rate statistics.

The measure is performed
comparing two binary streams and evaluating the B.E.R. as 
\f[
\widehat{BER}=\frac{n}{N}
\f] where \f$n\f$ is the number of bit in error
and \f$N\f$ the total number of bits compared. The
first stream can be delayed by a prefixed amount of symbols \f$d\f$ to
align the two compared sequence.


With the method SetParameters() the user can specify the delay \f$d\f$ of the first binary stream, 
and optionally the minimum number of bit \f$n_{\min}\f$ that should be counted in order to declare 
the estimate "reliable".

The method  Run()  delays the first binary stream by \f$d\f$ and  compares it with the second binary stream
maintaining a error counter \f$n\f$ and a bit counter \f$N\f$.

If the delay is specified to be a negative quantity  the BER meter starts in an "auto-align mode", In
auto-align mode the Run() method compares the second binary stream with all possible delayed version of the first
binary stream, with a delay \f$D\f$ ranging from 0 to \f$-d\f$. 

The method SetSoft() allows to pass the second binary stream as a sequence of soft values \f$\lambda\f$, as those
provided by soft decoders. 

Optionally, using the method SetFrameSize(), the user can specify the length of a frame \f$N_f\f$ 
so that the BER meter can perform also Frame Error Rate measures. 
The method requires to specify also the frame starting bit with respect to
the first binary stream. A frame is declared in error if at least one of its bit is in error.

The method IsReliable() returns a flag indicating if the BER meter has a reliable estimate of 
the bit (and frame if applicable) (\f$n\geq n_{\min}\f$).

For an example of its use see e.g. the test program "test_PCCC.cpp".

\author Guido Montorsi
*/
class BER_meter  
{
public:
	BER_meter();
	virtual ~BER_meter();

	//! Set the BER parameters
	void SetParameters(
		const int delay,	//!< Delay of second stream (if negative go to autolign mode)
		const int n=30		//!< Number of errors to declare realiable results
		);

	void SetSoft(bool s=true){soft=s;}	//!< Accepts soft information on second stream.

	//! Set a Frame size for FER statistics and specifies the relative offset 
	void SetFrameSize(
		const int fsizein,	//!< Frame size.
		const int foffin=0	//!< Starting offset with respect to the first stream.
		)
	{fsize=fsizein;foff=foffin;ferr=0;}



	//! Compare two bit streams. Return 1 if some error are present
	/*! The second stream can be represented as LLR if the SetSoft() method 
	has been called. */
	int Run(
		const int sizein,		//!< Number of bits. 
		const int* ref,			//!< First stream (source).
		const int* decoded		//!< Second bit stream (decoded).
		);



	bool IsReliable();	//!< Return true if enough bit or frame errors have been counted.
	void Reset();		//!< Reset the error counters

	//! Display the BER statistics on specified output.
	void Display(FILE* stream=stdout //!< Output stream
		); 

	//! Display the BER statistics on a single line of the specified output.
	void Display_on_Line(FILE* stream=stdout); 


	bool autoalign; //!< Find automatically the delay of second stream
	bool aligned;   //!< The BER meter is aligned and counting.
	int  nerrmin;	//!< Minimum number of errors to declare reliable measure.

	__int64 err;		//!< Number of bit errors
	__int64 nbit;		//!< Number of  bits in the statistic
	__int64 ferr;		//!< Number of frame errors
	__int64 nframe;		//!< Number of processed frames
	__int64 *stat;

	/** \cond INTERNAL */
	//! Carry on a statistic on the number of errors per frame (requires SetFrameSize).
	void StatErrors(
		const int minerr=1,	 //!< Minimum number of errors per frame considered in statistic.
		const int maxerr=10  //!< Maximum number of errors per frame considered in statistic.
		);

	//! Display the statistics on the number of errors per frame (see StatErrors()).
	void DisplayErrStat(
		FILE* stream=stdout,	//!< Output stream.
		bool singleline=false	//!< Flag to print the results on a single line
		);

	/** \endcond */

private:
	short* buffer; //!< Buffer for the delay
	int	*ber;		//!< Array of counters for autoalignment
	int	bufsize;	//!< size of delay buffer
	int	now;		//!< pointer to buffer
	int	delay;		//!< Delay between sequences
	__int64 bitcount;	//!< Number of bit from beginning of operation
	int minerr;
	int maxerr;

	bool soft;

	// FER measures
	int fsize;		//!< Size of frame. 
	int foff;		//!< Offset of starting frame
	int err_f;		//!< flag for erroneous frame
};
