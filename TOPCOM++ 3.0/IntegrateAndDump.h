// IntegrateAndDump.h: interface for the IntegrateAndDump class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class IntegrateAndDump
*/
/*! \ingroup DSP 
\brief Integrate and dump filter.

Through the SetParameters() method, the user specifies the
integral period \f$n\f$ expressed in number of samples. The output
signal is evaluated as:
\f[
s_{out}[i]=\sum_{\lfloor\frac{i}{n}-1\rfloor
n}^{\lfloor\frac{i}{n}\rfloor n}{s_{in}[i]}.
\f]

  For an example of its use see e.g. the test program "test_miscellanea.cpp".

\author Gabriella Bosco
*/
class IntegrateAndDump  
{
public:
	IntegrateAndDump();
	virtual ~IntegrateAndDump();

	//! Set teh main parameters of the block
	void SetParameters(
				const int int_period,		//!< Period, expressed in number of samples				
				const bool iscomplex=true	//!< Is the signal complex? (default=true)
				);
	//! Reset the integrator
	void Reset();

	//! Run the filter
	void Run (const int ntics,		//!< Number of input samples
			const double *Input,	//!< Input signal 
			double *Output			//!< Output signal
			);

#ifdef CTOPCOM
	void Run (const int ntics,		//!< Number of input samples
		const cmplx *Input,	//!< Input signal 
		cmplx *Output			//!< Output signal
		)
	{
			Run (ntics,		//!< Number of input samples
				(const double *)Input,	//!< Input signal 
				(double *)Output			//!< Output signal
				);
	}

#endif
private:
	bool iscmplx;		//!< Is the signal complex?
	int now;			//!< Time index
	int cont;			//!< Counts the samples in the integrator
	int period;			//!< Period, expressed in number of samples		
	double intr,		//!< Integral of the real part of the input signal
		inti;			//!< Integral of the imaginary part of the input signal
};

