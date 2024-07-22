// M_and_M_Timing.h: interface for the M_and_M_Timing class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
/*! \file 
\brief Declaration of class M_and_M_Timing
*/
/*! \ingroup Synchronization
\brief Closed-loop Muller & Mueller timing synchronizer.

This class implements a M_and_M_Timing closed loop timing recovery circuit. 
The sequence of samples \f$x_k\f$ at sample rate roughly \f$2/T\f$ enters an
interpolator AtoD whose sampling instant \f$\tau\f$ is controlled by the M_and_M_Timing loop. 
The output of the interpolator \f$y_n\f$ enters a shift register with three taps 
whose contents are sampled on even samples to generate the error signal.
The time synchronized output coincides with the even samples of \f$y_n\f$.

If not present outside, a match filter can be inserted on the sampled sequence to reduce the noise.

The method SetParameters of the class requires to specify the nominal value of the samples 
per period provided to the circuit (usually 2 but it can be more that that), 
the value of the updating step \f$\gamma\f$ and optionally the
pointer to the desired interpolator. 
If the interpolator is not specified the class creates a cubic interpolator with
oversampling factor 128.

Optionally the method SetFilter allows to insert a filter before entering the interpolator.

The Run method accepts the asynchronous sequence \f$x_k\f$ and provides the symbol synchronous 
sequence \f$y_{n/2}\f$. The
number of tics represents the number of processed input samples while the method returns the 
number of generated output
samples.

For an example of its use see e.g. the test program "test_Synchronization.cpp"
*/
#include "AtoD.h"
#include "Filter.h"
class M_and_M_Timing  
{
public:
	M_and_M_Timing();
	virtual ~M_and_M_Timing();

		//! Set the parameters of teh block
	void SetParameters(const double ns, //!< Nominal value of the samples  per period provided to the circuit 
		const double gamma,			//!< Value of the updating step 
		const double rho=0.,			//!< Value of the updating step 
		AtoD* AD=0					//!< Pointer to the desired interpolator
		);
	
	//! Insert a filter before entering the interpolator
	void SetFilter(Filter* MFin=0	//!< Pointer to a Filter
		);

	//! Run the synchronizer
	/*! The method return the number of generated output samples.*/
	int Run(const int tics,		//!< Number of processed input samples
		const double *inp,		//!< Input asynchronous sequence 
		double* out,			//!< Output synchronous sequence 
		const double* ref=0		//!< Output synchronous sequence 
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,		//!< Number of processed input samples
		const cmplx *inp,		//!< Input asynchronous sequence 
		cmplx* out,				//!< Output synchronous sequence 
		const cmplx* ref=0		//!< Output synchronous sequence 
		)

	{
		return (Run(tics,(const double*)inp,(double*) out,(double*) ref));
	}
	
#endif

	void Reset(){tau=0;}
	double tau;
	double e2;

private:
	void SetInterpolator(AtoD* AD=0);
	int nsamp; // Number of samples;
	double gamma;
	double rho;
	double *buff;
	double *newsamp;
	double *oldsamp;
	AtoD *Interp;
	Filter* MF;
	double N;	// Nominal sample per symbol at the input

};
