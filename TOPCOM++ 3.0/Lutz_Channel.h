// Lutz_Channel.h: interface for the Lutz_Channel class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "Gaussian_Process.h"

/*! \file 
\brief Declaration of the class  Lutz_Channel
*/
/*! \ingroup Channels
\brief The Lutz Channel 

	This class allows to generate a complex  process according to the Lutz model.
	In this model the channel can be in a bad state and in a good state 
	according to a transition diagram that is modeled through a two state Markov chain.
	
	The good state channel is a Rician channel with given rice factor
	and doppler frequency for the diffusive part.
	
	The bad state channel is a Rayleigh channel with average power modulated by a
	lognormal process. The parameters are
	the doppler frequency (the same as in the good channel) the variance of the lognormal
	process and the bandwidth of the lognormal process. 
	The spectrum of the lognormal process is a low-pass single pole spectrum, 
	while the spectrum of the fading is a classical Jake's spectrum.

For an example of its use see e.g. the test program "test_random_process.cpp".

\author Guido Montorsi

	*/
class Lutz_Channel  
{
public:
	Lutz_Channel();
	virtual ~Lutz_Channel();

	/*! Set the main parameters of the model */
	void SetParameters(
		const double dt,		//!< Time interval of generated process
		const double Rice,		//!< Rice factor for good channel [dB]
		const double fD,		//!< Doppler frequency 
		const double banddif,	//!< Bandwidth of lognormal shadowing  
		const double devdif,	//!< Standard deviation of lognormal shadowing 
		const double pgb,		//!< Probability of transition from good to bad channel
		const double pbg		//!< Probability of transition from bad to good channel
		);
	
	//! Generates complex samples of the Lutz process. 
	void	Run(const int tics,		//!< Number of generated samples
		double* proc				//!< Output samples of the Lutz process
		);
	
	//! Generates samples of the Lutz process and multiplies them for an input signal.
	void	RunProd(const int tics, //!< Number of generated samples.
		const double* input,		//!< Input process.
		double *output				//!< Output process.
		);

#ifdef CTOPCOM
	void	Run(
		const int tics,		//!< Number of generated samples
		cmplx* proc			//!< Output samples of the Lutz process
		)
	{
		Run(tics,(double*) proc);
	}
	void	RunProd(
		const int tics, //!< Number of generated samples.
		const cmplx* input,		//!< Input process.
		cmplx *output				//!< Output process.
		)
	{
		RunProd(tics,(const double*) input,(double *)output);
	}
#endif

private:
	double dt;			//!< Time interval of generated process
	double Rice;		//!< Rice factor (Ratio between LOS and diffusive power)
	double pgb;			//!< Bandwidth of LOS lognormal shadowing 
	double pbg;			//!< Standard deviation of  LOS lognormal shadowing 
	double fD;
	double banddif;		//!< Bandwidth of diffusive lognormal shadowing 
	double devdif;		//!< Standard deviation of diffusive lognormal shadowing 
	Gaussian_Process *diffusive;
	Gaussian_Process *shaddif;
	
	double *a,*b;
	double A1,A2,gain;
	int seed;
	int state;
	
	//	Power_Meter PM[3];
	
};
