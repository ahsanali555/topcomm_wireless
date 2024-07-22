// Equalizer.h: interface for the Equalizer class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdlib.h>
#include <stdio.h>
#include "Demodulator.h"
#include "pilot.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of class Equalizer */
/*! \ingroup Modems
\brief Adaptive Linear Equalizer.

  The block implements a general linear adaptive equalizer based on the stochastic 
  gradient algorithm.
  Minimum Mean square error (MMSE) equalizer and Decision Feedback equalizer (DFE)
  can be implemented with suitable initialization of the block.

  Fractionally spaced equalizer can be implemented with suitable initialization of the block.

  Initially, all taps coefficient are initialized to zero (Reset). 

  For an example of its use see e.g. the test program "test_equalizer.cpp".

\author Guido Montorsi
*/

class Equalizer  
{
public:
	Equalizer();
	virtual ~Equalizer();
	//! Set the main parameters 
	void SetParameters(
		Modulator *mod,			//!< Reference Memoryless Modulator
		const int Nf,			//!< Number of Feedforward symbols
		const double alphaf,	//!< Updating step FF
		const int nsin=1,		//!< Number of samples per symbol(>1 for fractionally spaced equalizers)
		const int Nb=0,			//!< Number of Feedback symbols  (>0 for DFE)
		const double alphab = 0,	//!< Updating step FB
		const double decw = 0.	//!< decision weight 
		);

	//! Set the length of the training sequence.
	/*! Returns the pointer to the training sequence.
	The training sequence is allocated and defined inside the block */
	double*  SetTraining(
		const int length,		//!< Length of the training sequence (number of symbols)
		int seed = 128901		//!< Seed for the generation of training sequence
		);

	//! Set the length of the training sequence.
	/*! Returns the pointer to the training sequence.
	The training sequence is allocated and defined inside the block */
	double*  SetTraining(
		const int length,		//!< Length of the frame
		const double* refpil	//!< Pilot sequence
		);

	//! Activates an optional PLL running after MMSE
	/*! Returns the pointer to the training sequence.*/
	void ActivatePLL(const double alphap, const double decwp = 1);

	//! Run the equalizer. 
	/*! The output buffer can coincide with the input buffer */
	void Run(
		const int ntics,		//!< Number of processed input samples
		const double* input,	//!< Input samples
		double* equal,			//!< Equalized samples (one per symbol)
		int* decoded=0,			//!< Decoded bits or symbols (optional)
		bool bin=true			//!< Decode bit or symbols?
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(
		const int ntics,		//!< Number of processed input samples
		const cmplx* input,	    //!< complex Input samples
		cmplx *equal,			//!< Equalized samples (one per symbol)
		int* decoded=0,			//!< Decoded bits or symbols (optional)
		bool bin=true			//!< Decode bit or symbols?
		)
	{
	Run(ntics,(const double*) input,(double*)equal,decoded,bin);
	return;
	}
#endif

	//! Print current taps coefficients 
	void Display(FILE *file=stdout);

	//! Reset taps and tap coefficients of equalizer.
	void Reset();
	double MSE;			//!< Current value of Mean square error, averaged on the last block of samples.
	double theta;       //!< Current value of phase estimate.
	int time;
	double power;

private:
	Modulator* mod;			//!< Reference modulator
	Demodulator* demod;
	int ns;				// Number of samples per period (fractional equalizer)
	double* f;          // Coefficients of filter
	double* b;          // Coefficients of filter
	double* linef;	    // Tapped delay forward line
	double* lineb;	    // Tapped delay backward line
	int nowf;			// Pointer to forward delay line
	int nowb;			// Pointer to backward  delay line
	int Nf;				// Number of Taps in the delay line
	int Nb;				// Number of Taps in the delay line
	int Ntrain;			// Length of training sequence
	double alpha;       // Upating FF factor
	double alphab;      // Upating FB factor
	double* training;   // Training sequence
	double decw;		// Decision weight for updating step
	double alphap;		// updating step of PLL
	double decwp;		// Decision weight for updating step in PLL
};

