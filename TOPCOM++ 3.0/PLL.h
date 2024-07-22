// PLL.h: interface for the PLL class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdlib.h>
#include <stdio.h>
#include "Demodulator.h"
#include "Pilot.h"
#include "Delay.h"
/*! \file 
\brief Declaration of class PLL */
/*! \ingroup Modems Synchronization
\brief Adaptive Linear PLL.

  The block implements a general PLL for phase recovery.

 

\author Guido Montorsi
*/

class PLL  
{
public:
	PLL();
	virtual ~PLL();
	//! Set the main parameters 
	void SetParameters(
		const Pilot* pilot,		//!< reference Pilot
		double alpha,			//!< updating step
		const Modulator* mod,	//!< Reference modulator for decisions
		const double wd			//!< Weight of data decision
		);


	//! Run the equalizer. Return the number of equalized data 
	/*! The output buffer cannot coincide with the input buffer */
	int Run(
		const int ntics,		//!< Number of processed input samples
		const double* inp,		//!< Input samples
		double* out,			//!< Equalized  output samples 
		double* phase=0,		//!< Equalized  output samples 
		bool strip= false		//!< Strip the pilots from sequence?
		);

	//! Run the equalizer. Return the number of equalized data 
	/*! The output buffer cannot coincide with the input buffer */
	int RunInterp(
		const int ntics,		//!< Number of processed input samples
		const double* inp,		//!< Input samples
		double* out,			//!< Equalized  output samples 
		double* phase=0,		//!< Equalized  output samples 
		bool strip= false		//!< Strip the pilots from sequence?
		);


	void Reset();	
	void Display(FILE* file=stdout);

	double alpha;
	double wd;
	bool off;

private:
	const Pilot* refpil;
	const Modulator* refmod;			//!< Reference modulator
	Demodulator* demod;
	int time;
	double theta;
	int npp;		// Counter of pilots positions
	int npi;		// Counter of pilots inserted


	// For interp 
	double th;
	double oth;
	double dth;
	bool final;
	bool first;
	double  RR;
	double II;
	Delay* delay;
};

