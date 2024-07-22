// Channel_Estimator.h: interface for the Channel_Estimator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Pilot.h"
#include "Delay.h"
#include "Filter.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class Channel_Estimator
*/
/*! \ingroup Synchronization
\brief Channel estimator based on pilots symbols (see also Pilot).

 This class implements the
functionality of a data aided channel estimator. 
The channel is assumed to be flat in frequency so that the estimated
channel is a complex random process that we will denote with \f$\alpha(t)\f$. 

The structure of the pilot stream is passed to the block through a pointer to the correspondent 
Pilot block.

The estimation of the channel process \f$\alpha\f$ is performed through correlation with the known 
pilot sequence \f$p(t)\f$ (whose timing is also assumed to be known) and filtering.
The type of filtering can be specified by the user through the methods SetIntegrateAndDump that specifies that
filtering is performed by integrate and dump, SetSlidingWindow that specifies a sliding window filter,
SetOnePole which  uses a one pole IIR filter and finally with the method SetEstimationFilter that allows the
user to choose the Filter. The correlation is performed only in the positions where the pilots are inserted, the
remaining positions in the stream are filled with the last available estimation.

For an example of its use see e.g. the test program "test_channel_estimation.cpp".

\author Guido Montorsi
*/
class Channel_Estimator  
{
public:
	//! Pass the  structure of the pilot stream
	void SetParameters(const Pilot* a		//!< Pointer to the reference Pilot block
		);

	//! Set the estimation filter type to Integrate And Dump
	void SetIntegrateAndDump(const int wind	//!< Time duration of the filter (in number of samples)
		);

	//! Set the estimation filter type to Sliding Window
	void SetSlidingWindow(const int wind	//!< Width of the window (in number of samples)
		);

	//! Set the filter type to one-pole IIR filter
	void SetOnePole(const double beta,		//!<
		const double alpha=0
		);

	//! Set the filter type to a user-defined one
	void SetEstimationFilter(Filter * filin,	//!< Pointer to a Filter	
		const int delay							//!< Delay of filter
		);

	//! Estimates the channel from noisy received samples
	int Run(const int ntics,					//!< Number of input samples
		const double* input,					//!< Input signal
		double* alpha,							//!< Estimated channel values
		double *output=0						//!< Equalized output (optional)
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,		//!< Number of input samples.
		const cmplx* input,			//!< Input complex signal
		cmplx*	alpha,				//!< Estimated channel value
		cmplx* output=0			//!< Output complex signal	
		)
	{
		Run(ntics,(const double*)input,(double*)alpha,(double*)output);
		return;
	}
#endif
	
	//! Reset the internal timing of the block
	void Reset(int off=0				//!< Initial offset of internal clock
		);

	Channel_Estimator();
	virtual ~Channel_Estimator();

private:
	int time;		// Internal clock
	int npp;		// Counter of pilots positions
	int npi;		// Counter of pilots inserted
	const Pilot* refpil;	// Pilot
	double corr[2];	// Last corrrelation result
	double a[2];	// Temporary

	int delay;		// delay of estimate process
	Filter* filter;	// Estimation filter
	int filtype;    // Type of filtering
	int nsamples;
	double window;

	double gain,beta;

	/* For removing overlapping pilots */
	Delay *deldat;  // Delay line for RX signal
	Delay *delpil;	// Delay line for pilot;
	Delay *delcor;	// Delay line for corr


};

