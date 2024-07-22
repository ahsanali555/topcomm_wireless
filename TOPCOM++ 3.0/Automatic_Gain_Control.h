// Automatic_Gain_Control.h: interface for the Automatic_Gain_Control class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include "Filter.h"
#include "ctopcom.h"

/*! \file
\brief Declaration of the Automatic_Gain_Control class 

*/

/*! \ingroup Synchronization
   \brief Performs automatic gain control, equalizing the sgnal power to a reference level.

The input signal samples \f$x_n\f$ , which can be real or complex, are multiplied by a constant \f$g\f$ whose value is
controlled by a feedback loop. To generate \f$g\f$ the power of samples \f$|x_n|^2\f$ is computed and compared with the
reference value obtaining

\f[ \Delta_n =P_{\rm ref}-|x_n|^2\f]

The gain \f$g\f$ is then updated according to the recursion
\f[
g^{(n+1)}=g^{(n)}+\gamma\Delta_n,
\f]
which represents a simple one pole filter (accumulator). An optional additional filter can be specified by the user if desired.


The method SetParameters() is used to specify the reference power level, the variable \f$\gamma\f$ and the type of input
signal (complex or real).

The optional filter can be added with the method SetFilter(). The method Freeze() can be used at any time to freeze
the current gain of the AGC or to un-freeze it.


For an example of its use see e.g. the test program "test_phase_jitters.cpp".

\author Guido Montorsi
*/  
class Automatic_Gain_Control  
{
public:
	//! Set the main parameters of block 
	void SetParameters(
		const double reflev, //!< Reference power level.
		const double gamma,	 //!< Updating step of loop (number of samples).
		bool complex=true	 //!< Flag for complex input.
		);
	//! Add an optional filter to the control block
	void SetFilter(Filter* filtin //!< Filter used to smooth the error signal
		);

	//! Run the block
	/*! The output buffer can coincide with the input buffer */
	void Run(const int tics, //!< Number of input samples.
		const double *inp,   //!< Input signal.
		double* out			 //!< Output signal
		);

	//! Display status
	void Display(FILE* file = stdout) const;


#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,		//!< Number of input samples.
		const cmplx* inp,			//!< Input complex signal
		cmplx* out					//!< Output complex signal	
		)
	{
		Run(tics,(const double*)inp,(double*) out);
		return;
	}		

#endif

	//! Freeze the AGC to current gain (or unfreeze)
	void Freeze(const bool fr=true	//!< Flag for freezing (or unfreezing)
		){freeze=fr;}

	Automatic_Gain_Control();
	virtual ~Automatic_Gain_Control();
	double gain;			//!< Current value of gain
private:
	Filter* filt;			//!< Optional filter.
	double reflev;			//!< Reference amplitude level.
	bool complex;			//!< Flag for complex input.
	double gamma;			//!< Updating step of loop.
	bool freeze;
	double mean;
	double stdev;
};
