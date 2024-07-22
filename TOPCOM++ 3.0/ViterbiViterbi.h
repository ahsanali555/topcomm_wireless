// ViterbiViterbi.h: interface for the ViterbiViterbi class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Modulator.h"
#include "Filter.h"
#include "Pilot.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class ViterbiViterbi
*/
/*! \ingroup Synchronization 
\brief  Viterbi & Viterbi algorithm for phase recovery

  The class implements both feedback and feedforward algorithms.
  Parameters are the (linear) Modulator used at the transmitter side, 
  the power used for the non-linear transformation and the length 
  of the observation window.

For an example of its use see e.g. the test program "test_phase_synchronization.cpp".

\author Guido Montorsi

*/


class ViterbiViterbi  
{
public:
	ViterbiViterbi();
	virtual ~ViterbiViterbi();

	//! Set the main parameters for the Viterebi&Viterbi algorithm
	void SetParameters(
		const int L,	//!< [in] Length of observation window
		const int m,	//!< [in] Phase ambiguity (fraction of 2pi) of modulation
		const int pow=4	//!< [in] Power used for the non linearity
		);

	//! Set the loop filter of the V&V with feedback
	void SetFeedback(
		const double A,	//!< [in] Gain of error signal
		const int fil=1 //!< [in] Order of filter used in the loop (default to one pole)
		);

	//! Run the V&V algorithm. 
	/*! The phase estimate is accessible through the variable theta */
	void Run(const int tics,	//!< [in] Number of processed complex samples
		const double*input ,	//!< [in] Input samples
		double* out=0			//!< [out] Sequence of phase compensated samples [optional]
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	//! Run the V&V algorithm. 
	/*! The phase estimate is accessible through the variable theta */
	void Run(const int tics,	//!< [in] Number of processed complex samples
		const cmplx*input ,		//!< [in] Complex Input samples
		cmplx* out=0			//!< [out] Complex Sequence of phase compensated samples [optional]
		)
	{
	Run(tics,(const double*)input ,(double*)out);
	return; 
	}
#endif
	
	//! Reset the estimator
	/*! The counter for observed samples and the accumulated values for
	real and imaginary parts are reset to zero. The observation window is then
	aligned to the next incoming sample */ 
	void Reset();

	double theta;		//!< Current phase estimate

private:
	const Pilot* refpil;	
	Filter* filt;
	bool feedback;
	double A;
	double damp;
	int pow;
	int L;
	int time;
	double X,Y,ccos,ssin;
	int m;

};

