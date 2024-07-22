// Histogram.h: interface for the Histogram class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "Include.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of the class Histogram
*/
/*! \ingroup Measurements 
\author Gabriella Bosco, Guido Montorsi
\brief  Store the histogram of a signal.

  The Histogram class implements a block generating the histogram of a signal 
  in the range \f$\{h_{min}, h_{max}\}\f$.
  If \f$s_{in}[i]\f$ is the input signal, the following samples are used to evalute the histogram: 
\f$ s_{in}[n_{offset}+k\cdot \Delta n_s]\;\;\;\;k=0,1,...\f$, 
where \f$\Delta n_s\f$ is the sampling period expressend in number of samples.
The number of samples which fall out from the range \f$\{h_{min}, h_{max}\}\f$ is counted 
by the two integer variables
\f$h_{underflow}\f$ and \f$h_{overflow}\f$ and stored in the output file as well (in first and last position).

The user specifies:
- The name of output file in which the histogram is stored.
- The minimum value of histogram range \f$h_{min}\f$.
- The maximum value of histogram range \f$h_{max}\f$.
- The number of bins.
- The offset \f$n_{offset}\f$ (expressed in number of samples).
- The sampling period \f$\Delta n_s\f$ (expressed in number of samples).

  For an example of its use see e.g. the test program "test_measurements.cpp".

*/

class Histogram  
{
public:
	Histogram();
	virtual ~Histogram();

	//! Set the main parameters of the block
	void SetParameters(
				const char* name,	//!< Name of output file
				const double min,	//!< Minimum value of histogram range
				const double max,	//!< Maximum value of histogram range
				const int nbin=100,	//!< Number of bins (default=100)
				const int offset=0,	//!< Offset, expressed in number of samples (default=0)	
				const int sampling_period=1  //!< Sampling period of data (default=1)	
				);

	//! Reset the histogram storage
	void Reset();

	//! Start the storage of the signal samples
	void Start(){start=true;nsamples=0;};

	//! Terminate the storage of the signal samples
	void Stop();

	//! Return the mean value of the stored samples
	double Mean(){return mean/nsamples;}

	//! Return the variance of the stored samples
	double Var() {return m2/nsamples - Mean() * Mean();}

	//! Return the abscissa yielding the specified value of cumulative
	double AbscissaAt(const double p//!< value of cumulative
		);
	//! Print the histogram in an output file
	void Print(FILE* file=0,		//!< Name of output file (default: "Histogram#.txt")
				bool line=false		//!< Output in a line.
				);
	
	//! Print the header of the output file
	void PrintHeader(FILE *f);
	//! Run the histogram block
    void Run(const int ntics,   //!< Number of input samples
			 double* Input		//!< Input signal or quantity
			 );
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	   void Run(const int ntics,   //!< Number of input samples
			 cmplx* Input		//!< Input signal or quantity
			 )
	   {
	   Run(ntics,(double*) Input);
	   }
#endif


	__int64 *occurrences;		//!< Number of occurrences of each bin
	__int64 overflows;			//!< Number of overflows (signal samples higher than \f$h_{max}\f$) 
	__int64 underflows;			//!< Number of underflows (signal samples lower than \f$h_{min}\f$) 
	__int64 nsamples;			//!< Total number of stored samples
	double mean;				//!< Mean value of stored samples
	double m2;					//!< Mean value of squared stored samples

private:

	bool start;
	double delta;
	double min;
	double max;

	int now;
	int step;
	int dim;
	FILE* file;
};

