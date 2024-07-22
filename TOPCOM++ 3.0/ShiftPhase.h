// ShiftPhase.h: interface for the ShiftPhase class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class ShiftPhase
*/
/*! \ingroup DSP 
\brief Offsets the phase of a complex discrete signal.

Through the SetParameters() method, the user specifies the phase shift
$\Delta\phi$ expressed in radians. If \f$s_{in}[i]\f$ is the complex
envelope of the input signal, the output signal \f$s_{out}[i]
\f$ is given by:
\f[
\begin{array}{l}
Re\{s_{out}[i]\}=\cos(\Delta\phi)\cdot Re\{s_{in}[i]\}-\sin(\Delta\phi)\cdot Im\{s_{in}[i]\}\\
Im\{s_{out}[i]\}=\cos(\Delta\phi)\cdot Im\{s_{in}[i]\}+\sin(\Delta\phi)\cdot Re\{s_{in}[i]\}
\end{array}
\f]

For an example of its use see e.g. the test program "main_miscellanea.cpp".

\author Gabriella Bosco
*/
#include "ComplexGain.h"

class ShiftPhase  
{
public:
	//! Set the phase shift
	ShiftPhase(const double phase    //!< Phase shift in radians
		);
	virtual ~ShiftPhase();
	//! Run the pahse shifter
	void Run(const int ntics,		//!< Number of samples
			 const double* Input,	//!< Input signal
			 double* Output			//!< Output signal
			 );
#ifdef CTOPCOM
	void Run(const int ntics,		//!< Number of samples
		const cmplx* Input,	//!< Input signal
		cmplx* Output			//!< Output signal
		)
	{
		Run(ntics,		//!< Number of samples
			(const double*) Input,	//!< Input signal
			(double*) Output			//!< Output signal
			);
	}
#endif
private:
	ComplexGain* gain;				//!< Phase shift expressed as a complex number
};

