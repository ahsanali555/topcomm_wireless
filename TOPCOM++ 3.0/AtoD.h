// AtoD.h: interface for the AtoD class.
//
//////////////////////////////////////////////////////////////////////
/*! \file
\brief Declaration of the class AtoD and its interface functions 
*/
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \ingroup DSP
\brief Generic interpolator and resampling block.


This blocks implements all the functions typical of a real Analog to Digital converter.
Upsampling, downsampling , re-sampling and sampling controlled by a timing signal affected
by an arbitrary error can be simulated through this block.
The block manage both complex and real signals.

See the class Sampler for a simpler decimator block.

Interpolation is performed through an interpolating filter that can be specified by the user 
through the constructor of the class.
The interface methods Poly_AtoD(), Linear_AtoD() and Raised_Cosine_AtoD()
allows to construct A/D converters with typical interpolators.


To speedup the interpolation process, the input signal is oversampled by a factor specified 
by the user and then resampled to the desired sampling instants (up to the oversampling precision).
The oversampling factor does not increase the complexity of the interpolating 
algorithm but only the memory requirements of the interpolating coefficients.

The sampling frequency \f$f_c\f$, a real number specified in term of
number of output samples per input sample, can be changed at any
time through SetSamplingFrequency(). Alternatively the user can
specify the instantaneous frequency signal through a vector that is
passed together with the input signal and with the same rate of it
to the Run() method. This allows to deal also with rapidly
changing timing signals. The timing offset \f$\Delta\f$ is the amount of
input samples (not necessarily integer) that are skipped before the
first output sample. At the beginning \f$\Delta\f$ is set to 0 and the
user can change the timing offset using ShiftSamplingTime(). 
A positive value shift means that the next
output will be delayed by the specified amount, while a negative
value means that the next output will be anticipated by the
specified amount of input samples.

As a result of the time varying nature of the AtoD timing signal,
the AtoD  block computes a number of output samples that  can change
from block to block. The number of computed output samples is
returned from the Run() method.

For an example of its use see e.g. the test program "test_AtoD.cpp".

\author Guido Montorsi
*/
class AtoD  
{
public:

	//! Create an AtoD. 
	/*!
	 The object is defined through the interpolating filter coefficients and the oversampling factor.
	*/
	AtoD(const int oversampling,  //!< Oversampling of the input signal
		 const int ntaps,		  //!< Number of taps of interpolating filter
		 const double * interp,	  //!< Coefficients  of interpolating filter (their number is oversampling x ntaps)
		 bool iscmplx=false		  //!< Flag to indicate if complex or real AtoD
		 );	  


	virtual ~AtoD();

	//! Samples the input signal, returns the number of generated samples 
	int	Run(
		const int nsamples,			//!< Number of processed input samples
		const double* input,		//!< Input samples
		double* output,				//!< output samples.
		double* freq=0				//!< Optional instantaneous sampling frequencies for all input samples
		);	

	#ifdef CTOPCOM
	int	Run(
		const int nsamples,			//!< Number of processed input samples
		const cmplx* input,		//!< Input samples
		cmplx* output,				//!< output samples.
		double* freq=0				//!< Optional instantaneous sampling frequencies for all input samples
		)
	{
		return Run(nsamples, (const double*)input, (double*)output, freq);
	}
	#endif
	

	//! Shift the sampling time with respect  to the current one. 
	void    ShiftSamplingTime(
		const double shift			//!< Shift of the sampling time
		);

	//! Set the current sampling frequency (normalized to input sampling rate). 
	void    SetSamplingFrequency(
		const double frequency		//!< Sampling frequency  (output samples per input sample)
		){fr_ist=frequency;};

	//! Set the current sampling frequency and time offset. 
	void    SetSamplingFrequency(const double frequency, //!< Sampling frequency  (output samples per input sample)
			const double dphase						   	//!< Shift of the sampling time
		);

	//! Return the delay introduced by the A/D block
	double	GetDelay()
	{
		return (double)(ntaps/2)+phist/fr_ist;
	}



	//! Reset the block 
	void	Reset();
	double  phist;			//!< Instantaneous phase
	double fr_ist;			//!< Istantanuous frequency of sampler
	
private:
	int oversampling;		//!< Oversampling of interpolating filter
	int ntaps;				//!< Number of taps interpolating filter
	double* interp;			//!< Coefficients of the interpolating filter
	double* filter;			//!< Filter tapped delay line
	int pointer;			//!< pointer for TDL 
	bool iscmplx;			//!< The interpolator works on complex numbers?

friend	AtoD* Poly_AtoD(const int oversampling,const int degree,bool iscmplx);
friend AtoD* Linear_AtoD(const int oversampling,bool iscmplx);
friend AtoD* Raised_Cosine_AtoD(const int oversampling,const double band,double rolloff,bool iscmplx);

};


	/*! \ingroup Interface
\brief Generate an A/D converter with a polynomial interpolator 
*/

AtoD* Poly_AtoD(const int oversampling,		//!< Oversampling of the input signal
				const int degree=3,			//!< Degree of the polynomial (default to cubic interpol)
				bool iscmplx=false			//!< Is the signal complex?	 
				);

/*! \ingroup Interface
\brief Generate an A/D converter with a linear interpolator 
*/

AtoD* Linear_AtoD(const int oversampling,	//!< Oversampling of the input signal
				  bool iscmplx=false		//!< Is the signal complex?	 
				  );


/*! \ingroup Interface
\brief Generate an A/D converter with a raised cosine interpolator 
*/


AtoD* Raised_Cosine_AtoD(const int oversampling, //!< Oversampling of the input signal
						 const double band=0.8,	 //!< Flat band of RC filter (<1)
						 double rolloff=-.1,	 //!< Roll-off of RC filter 
						 bool iscmplx=false		//!< Is the signal complex?	 
						 );
