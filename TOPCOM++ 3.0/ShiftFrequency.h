// ShiftFrequency.h: interface for the ShiftFrequency class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class ShiftFrequency
*/
/*! \ingroup DSP 
\brief Offsets the central frequency of a complex discrete signal.

The user specifies the product between the frequency shift \f$\Delta f\f$
and the sampling interval \f$\delta t_s\f$,
expressed as the ratio between two integers {\bf num} and {\bf den}.
The rational form of the normalized frequency allows
to implement the frequency shift with fast look-up tables with size equal to the denominator.


For an example of its use see e.g. the test program "main_highratetelemetry.cpp".

\author Guido Montorsi
*/
class ShiftFrequency  
{
public:
	ShiftFrequency();

	//! Set the frequncy shift
	ShiftFrequency(int num,			//!< Numerator of the product between the frequency shift and the sampling time
					int den,		//!< Denominator of the product between the frequency shift and the sampling time
					int offset=0	//!< Initial phase offset
						);
	virtual ~ShiftFrequency();

	//! Run the frequency shifter
	void Run(const int tics,	//!< Number of samples
		const double *input,	//!< Input signal
		double *output		//!< Output signal
		);

	void ChangeFrequency(const double f);
	void ChangePhase(const double phase);

	unsigned int num;	 //!< Numerator of the product between the frequency shift and the sampling time
	unsigned int den;	//!< Denominator of the product between the frequency shift and the sampling time
	unsigned int tempo;	//!< Time index
#ifdef CTOPCOM
	void Run(const int tics,	//!< Number of samples
		const cmplx *input,	//!< Input signal
		cmplx *output		//!< Output signal
		)
	{
		Run(tics,	//!< Number of samples
		(const double *)input,	//!< Input signal
		(double *)output		//!< Output signal
		);
	}
#endif
private:
	double* s;	//!< Vector containing the values of \f$\sin(2\pi\cdot i/{\rm den})\;\;\;i=1,..,{\rm den}\f$
	double* c;	//!< Vector containing the values of \f$\cos(2\pi\cdot i/{\rm den})\;\;\;i=1,..,{\rm den}\f$
};
