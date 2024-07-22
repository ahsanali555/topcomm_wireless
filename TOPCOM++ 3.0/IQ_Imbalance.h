// IQ_Imbalance.h: interface for the IQ_Imbalance class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class IQ_Imbalance
*/
/** \ingroup Channels
\brief   Generates an I/Q imbalance in phase, offset and amplitude.

The class implements the
impairments on the received signal due to the non perfectly balanced
amplitudes and  phases of the two sinusoids for the baseband
conversion of the RF signal. The effects are described through the
amplitude imbalance  \f$A\f$ between In-phase and Quadrature components
(expressed in dB), the relative phase imbalance \f$\varphi\f$, and the
complex offset \f$O=O_I+j O_Q\f$.

The two sinusoids are then assumed to be
   \f[
a_I(t)=2\sqrt{1+g}\cos(2\pi f_0t-\varphi/2)\f]
\f[a_Q(t)=2\sqrt{1-g}\sin(2\pi f_0t+\varphi/2)\f]
such that
\f[
\left.{\frac{1+g}{1-g}}\right|_{dB}=A \rightarrow g=
\frac{a-1}{a+1}\;\;a=10^{A/10}
\f]

The  I and Q signals are then obtained as follows:
  \f[I(t) =O_I + \sqrt{1+g}(x_I(t)\cos\varphi/2 + x_Q(t)\sin\varphi/2) \f]
  \f[Q(t) =O_Q +  \sqrt{1-g}(x_Q(t)\cos\varphi/2 + x_I(t)\sin\varphi/2)\f]

The method SetParameters()  can be used to set or modify \f$A,
\varphi, O_I\f$ and \f$O_Q\f$.

For an example of its use see e.g. the test program "test_miscellanea.cpp".

\author Guido Montorsi
*/
class IQ_Imbalance  
{
public:
	IQ_Imbalance();
	virtual ~IQ_Imbalance();	//! Simulation method
	void Run(const int tics,	//!< Number of processed input samples
		const double* inp,		//!< Input complex samples 
		double* out				//!< Output samples affected by I/Q imbalance
		);

	//! Set the main parameters
	void SetParameters(const double amp, //!< Amplitude imbalance [dB]
		const double phase=0,			//!< Phase imbalance (degrees)
		const double Ioff=0,			//!< Inphase offset
		const double Qoff=0				//!< Quadrature offset
		);

#ifdef CTOPCOM
	void Run(const int tics,	//!< Number of processed input samples
		const cmplx* inp,		//!< Input complex samples 
		cmplx* out				//!< Output samples affected by I/Q imbalance
		)
	{
		Run(tics,(const double*) inp,(double*) out);
	}
#endif
private:
	double ampI;
	double ampQ;
	double phase;
	double Ioff;
	double Qoff;
	double c,s;

};

