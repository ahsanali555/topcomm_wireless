// Sampler.h: interface for the Sampler class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class Sampler

*/
/*! \ingroup DSP
\brief Decimates a discrete sequence

Decimates the input signal \f$s_{in}[i]\f$ with a sampling period \f$n\f$ and an 
offset \f$k\f$ (both
expressed in number of samples), defined by the user through the SetParameters() method. The output signal is given
by:
\f[
s_{out}[i]=s_{in}\left[k+i\cdot n\right].
\f]

See AtoD for more complex interpolating functions

For an example of its use see e.g. the test program "main_highratetelemetry.cpp".

\author Guido Montorsi
*/
class Sampler  
{
public:
	Sampler();
	virtual ~Sampler();

	//! Set the main parameters of the block	
	void SetParameters(const int nsin,				//!< Period, expressed in number of samples	
					const int offsetin,				//!< Offset, expressed in number of samples	
					const bool icmplx=false			//!< Is th einput signal complex? (default:false)
					);
	//! Run the decimator, return the number of samples.
	int Run(const int tics,		//!< Number of input samples
			const double* input ,	//!< Input signal 
				double *output		//!< Output signal
				);
#ifdef CTOPCOM
	int Run(const int tics,		//!< Number of input samples
		const cmplx* input ,	//!< Input signal 
		cmplx *output		//!< Output signal
		)
	{
		return Run( tics,		//!< Number of input samples
			(const double*) input ,	//!< Input signal 
			(double *)output		//!< Output signal
			);
	}
#endif
private:
	int ns;			//!< Period, expressed in number of samples	
	int offset;		//!< Offset, expressed in number of samples	
	int time;		//!< Time index
	int icmplx;		//!< Is the signal complex?
};
