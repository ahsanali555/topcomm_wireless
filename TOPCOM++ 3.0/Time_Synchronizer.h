// Time_Synchronizer.h: interface for the Time_Synchronizer class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "AtoD.h"
#include "Delay.h"
#include "Filter.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class Time_Synchronizer
*/
/*! \ingroup Synchronization
\brief Oerder & Meyr open loop timing synchronizer with cubic interpolator.

The method SetParameters() of the class requires to specify the nominal value of the samples 
per period provided to the algorithm (usually 4 but it can be more than that), 
the value of the updating observation window \f$L\f$ and optionally
the pointer to the desired interpolator. 
If the interpolator is not specified the class creates a cubic interpolator
with oversampling factor 128.

For an example of its use see e.g. the test program "test_synchronization.cpp" 

\author Guido Montorsi
*/
class Time_Synchronizer  
{
public:
	Time_Synchronizer();
	//! Set the parameters of the synchronizer
	void SetParameters(const int N,	//!< Nominal value of the samples per period provided to the algorithm
		const int L0,				//!< Value of the updating observation window
		AtoD* AD=0					//!< Pointer to the desired interpolator
		);
	//! Set the interpolator 
	void SetInterpolator(AtoD* AD=0 //!< Pointer to an AtoD used for interpolation
		);
	//! Set the filter
	void SetFilter(Filter* MFin=0	//!< Pointer to a Filter
		){MF=MFin;};

	void Display(FILE*file = stdout) const;

	void Reset() { sumI = sumR = 0.; time = 0; };
	virtual ~Time_Synchronizer();

	//! Run the synchronizer
	/*! The method return the number of generated output samples.*/
	int Run(const int tics, //!< Number of processed input samples
		const double *inp,	//!< Input asynchronous sequence 
		double* out			//!< Output synchronous sequence 
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	//! Run the synchronizer
	/*! The method return the number of generated output samples.*/
	int Run(const int tics, //!< Number of processed input samples
		const cmplx *inp,	//!< Complex Input asynchronous sequence 
		cmplx* out			//!< Complex Output synchronous sequence 
		)
	{
	return (Run(tics,(const double*)inp,(double*) out));
	}
#endif

	double tau;				//!< estimated offset
	AtoD* Interp;
	bool frozen;

private:
	int N; // Number of samples per symbol
	int L0; // Observation window
	int time;
	double sumR,sumI;
	double *c;
	double *s;
	Delay	*delay;
	Filter* MF;
	double mean;
	double stdev;
};
