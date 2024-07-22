// Gaussian_Process.h: interface for the Gaussian_Process class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "AWGN_Channel.h"
#include "Filter.h"
#include "AtoD.h"
#include "Spectrum.h"
/*! \file 
\brief Declaration of class Gaussian_Process
*/
/*! \ingroup Random_Sources 
\brief Generates a generic Gaussian process with user specified power spectral density.

  The block is realized
through the cascade of a white Gaussian generator with unitary power, followed by one 
or more parallel filters and by a
linear interpolator. The interpolation is used to efficiently generate processes that
have a bandwidth considerably smaller than the simulation bandwidth as is often the case 
(e.g. noise process or fading
processes).

Through the method SetParameters() the user specifies the  filters that are used to shape
 the spectrum of the
generated process and the upsampling ratio of the interpolator \f$f_u>1\f$.


The Run() method generates a sequence of Gaussian samples with a sampling rate that is 
\f$\frac{1}{f_u}\f$ times the
desired sample frequency. The generated samples are filtered and then linearly interpolated 
to the required sampling
frequency. The RunProd() uses the generated samples to modulate (multiply) an input signal provided by the user.

If \f$f_u\gg 1\f$ the complexity of the generation and filtering of noise samples becomes negligible with respect to the
overall complexity, which is then dominated by the interpolator.

In order to provide reliable samples with linear interpolation the bandwidth of the filters must be considerably
smaller than the Nyquist bandwidth.

The method NormalizePower()  normalizes the power of the generated
process to a desired value, considering also the gain introduced by
the filter and the interpolator.

If the spectrum of the desired process has components that cover a large 
frequency interval and with very different
amplitudes we suggest to represent the process as the sum of
 processes with different associated powers.

 For an example of its use see e.g. the test program "test_channel_estimation.cpp".

\author Guido Montorsi 
*/
class Gaussian_Process
{
public:
	Gaussian_Process();
	virtual ~Gaussian_Process();

	//! Set the parameters of the Gaussian process.
	void SetParameters(
		Filter *filter,				//!< Filter(s) to shape the autocorrelation function
		const double fratio=1.,		//!< Linear upsampling ratio
		const int nfilters=1		//!< Number of parallel filters
		);

	
	//! Scale the amplitude of process
	void SetGain(
		const double g //!< Factor for which the current noise power is multiplied
		)
	{
		gain*=g;
		if(iscmplx)whitenoise->Set_EsN0(1. / gain);
		else	   whitenoise->Set_EsN0(0.5/ gain);
	}

	//! Normalize the power of the process 
	void NormalizePower(
		const double pow=1.	//!< Power of the process
		);


	//! Reset the generator 
	void Reset();
	
	//! Generates samples of a Gaussian process. 
	void	Run(const int tics,		//!< Number of generated samples
			double* proc			//!< Output samples of the Gaussian process
			);

	//! Generates samples of a Gaussian process and multiplies them for an input signal.
	double*	RunProd(const int tics, //!< Number of generated samples.
		const double* input,		//!< Input signal.
		double *output				//!< Output process.
		);

#ifdef CTOPCOM
	//! Generates complex samples of a Gaussian process. 	
	/*! dettagli*/
	void Run(const int tics,		//!< Number of generated samples
			cmplx* proc				//!< Output samples of the Gaussian process
			)
	{
		Run(tics, (double*)proc);
		return;
	}

	cmplx*	RunProd(const int tics, //!< Number of generated samples.
		const cmplx* input,		//!< Input signal.
		cmplx *output				//!< Output process.
		)
	{
		return (cmplx*)RunProd(tics, (const double*) input,	(double *)output);
	}
#endif



	//! Measure the PSD of the generated process.
	void MeasurePSD(
		Spectrum* a,	//!< pointer to PSD measure block 
		bool meas=true	//!< measure PSD or not?
		);

	
	friend 
		Gaussian_Process* Jakes_Process(const double fd,const double Rs );
		
	friend 
		Gaussian_Process* One_Pole_Process(const double band,const bool icmplx);


	bool iscmplx;					//!< Is the random process complex?
	AWGN_Channel *whitenoise;		//!< White Gaussian generator(s)
	Filter		 *filter;			//!< Filter(s) Shaping the spectrum
	int nfil;						//!< Number of filters. 

private:
	bool wrap;
	double gain;
	bool	meas;
	Spectrum *sp;
	double fratio;  //!< Ratio between sample rate before and after interpolation
	double *noise;	//!< Gaussian noise buffer
	double time;

	double *oldp;
	double *newp;
	double tmp[2];
};

/*! Returns a complex Gaussian process with Jake's spectrum. */	
Gaussian_Process* Jakes_Process(const double fd,//!< Doppler frequency
		const double Rs 								//!< Simulation sampling rate					
		);

//! Returns a real Gaussian process with a low pass one pole spectrum */
Gaussian_Process* One_Pole_Process(const double band,		//!< Bandwidth of filter.
											const bool icmplx=false //!< Flag for complex input signals
								);






