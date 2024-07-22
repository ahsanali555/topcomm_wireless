// Loo_Channel.h: interface for the Loo_Channel class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "Gaussian_Process.h"
#include "Power_Meter.h"
/*! \file 
\brief Declaration of the class Loo_Channel
*/

/*! \ingroup Channels
\brief The model for Loo and Corazza fading channels 

This class allows to generate a complex  process according to the Loo model.

In this model the process is obtained as the sum of a diffusive component 
and a Line Of Sight (LOS) component.
The ratio between the power of the two components is determined by the Rice factor.

The diffusive component is modeled as a Rayleigh channel with Jake's spectrum
with average power modulated by a lognormal process. 
  
The LOS component is modeled as a constant with average power modulated 
by a lognormal process. 

The parameters are
the doppler frequency ,the variance and the bandwidth of the two lognormal processes
that modulate the power of LOS and diffusive component.

The spectra of the lognormal processes are  low-pass single pole spectra, 
while the spectrum of the fading is a classical Jake's spectrum.
			
\see Lutz_Channel()

For an example of its use see e.g. the test program "test_random_process.cpp".

\author Guido Montorsi
*/

class Loo_Channel  
{
public:
	Loo_Channel();
	
	/*! Set the main parameters of the model */
	void SetParameters(
		const double dt,		//!< Time interval of generated process
		const double Rice,		//!< Rice factor [dB] (Ratio between LOS and diffusive power)
		const double bandLOS,	//!< Bandwidth of LOS lognormal shadowing 
		const double devLOS,	//!< Variance of LOS lognormal shadowing 
		const double fD,		//!< Doppler frequency of diffusive component
		const double banddif=0,	//!< Bandwidth of diffusive lognormal shadowing 
		const double devdif=0		//!< Variance of diffusive lognormal shadowing 
		);
	
	//! Generates samples of the Loo process. 
	void	Run(const int tics,		//!< Number of generated samples
		double* proc			//!< Output samples of the Gaussian process
		);
	
	//! Generates samples of the Loo process and multiplies them for an input signal.
	void	RunProd(const int tics, //!< Number of generated samples.
		const double* input,		//!< Input signal.
		double *output				//!< Output process.
		);
	//! Display Channel parameters
	void	Display(FILE* =stdout);
	
	virtual ~Loo_Channel();

#ifdef CTOPCOM
	void	Run(
		const int tics,		//!< Number of generated samples
		cmplx* proc			//!< Output samples of the Gaussian process
		)
	{
			Run(tics,(double*) proc);
	}

	//! Generates samples of the Loo process and multiplies them for an input signal.
	void	RunProd(const int tics, //!< Number of generated samples.
		const cmplx* input,		//!< Input signal.
		cmplx *output				//!< Output process.
		)
	{
		RunProd(tics,(const double*) input,(double *)output	);
	}
#endif
private: 
	double dt;		//!< Time interval of generated process
	double Rice;	//!< Rice factor (Ratio between LOS and diffusive power)
	double bandLOS;	//!< Bandwidth of LOS lognormal shadowing 
	double devLOS;	//!< Standard deviation of  LOS lognormal shadowing 
	double fD;
	double banddif;	//!< Bandwidth of diffusive lognormal shadowing 
	double devdif;	//!< Standard deviation of diffusive lognormal shadowing 
	
	Gaussian_Process *diffusive;
	Gaussian_Process *shaddif;
	Gaussian_Process *shadLOS;
	
	double *a,*b,*c;
	double A1,A2,gain;
	
	//	Power_Meter PM[3];
};
