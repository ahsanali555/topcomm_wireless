// SampleAndHold.h: interface for the SampleAndHold class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class SampleAndHold
*/
/*! \ingroup DSP 
\brief Sample and hold block.

Through the SetParameters() method, the user specifies the
sampling period \f$n\f$ and the offset \f$k\f$, both expressed in number of
samples. The output signal is evaluated as:
\f[
s_{out}[i]=s_{in}\left[k+\left\lfloor\frac{i}{n}\right\rfloor \cdot n\right].
\f]

For an example of its use see e.g. the test program "test_miscellanea.cpp".

\author Gabriella Bosco
*/


class SampleAndHold  
{
public:
	SampleAndHold();
	virtual ~SampleAndHold();
	//! Set the main parameters of the block
	void SetParameters(
				const int int_period,		//!< Period, expressed in number of samples				
				const int offset=0,			//!< Offset, expressed in number of samples	
				const bool iscomplex=true	//!< Is the input signal complex? (default: true)
				);
	//! Reset the sample and hold
	void Reset(){cont=0;};
	//! Run the block
	void Run (const int ntics,				//!< Number of input samples
		const double *Input,				//!< Input signal 
		double *Output					//!< Output signal
		);

#ifdef CTOPCOM
	void Run (const int ntics,				//!< Number of input samples
		const cmplx *Input,				//!< Input signal 
		cmplx *Output					//!< Output signal
		)
	{
		Run (ntics,				//!< Number of input samples
			(const double *)Input,				//!< Input signal 
			(double *)Output					//!< Output signal
			);
	}
#endif
private:
	bool iscmplx;			//!< Is the signal complex?
	int now;				//!< Time index
	int cont;				//!< Counts the samples in the sample and hold block
	int period;				//!< Period, expressed in number of samples
	double memr,			//!< Real part of the stored sample
		memi;				//!< Imaginary part of the stored sample


};
