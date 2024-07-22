// Fading_Channel.h: interface for the Fading_Channel class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdio.h>
#include "FFT.h"
#include "Power_Meter.h"
/*! \file 
\brief Declaration of the class Fading_Channel and its associated interfaces
*/

/*! \ingroup Channels
\brief The multipath fading channel.

Generic Single Input Single Output multipath fading channel. All taps coefficients
are Gaussian random processes with Jake's spectrum. The delay of
each tap is kept constant during the simulation. Change of
the delay associated to a tap can be made by using the method
ChangeFadingTaps(). To improve the efficiency of the block, fading
processes are always generated at a sampling rate equal to 8 time
the normalized doppler frequency. Oversampling with linear
interpolation is performed to achieve the required sampling rate.

All processes of tap coefficients have zero mean. If one need 
to generate Rice processes, use the more general class MIMO_Fading_Channel instead.

For an example of its use see e.g. the test program test_CompCodeComb.cpp
and its interface functions.

\see  
		- Multipath_COST207_TU6()
		- Multipath_UMTSCase0(const double dt, const double fd)
		- Multipath_UMTSCase1(const double dt, const double fD)
		- Multipath_UMTSCase2(const double dt, const double fD)
		- Multipath_UMTSCase3(const double dt, const double fD)
		- Multipath_UMTSCase4(const double dt, const double fD)
		- Multipath_UMTSCase5(const double dt, const double fD)
		- Multipath_UMTSCase6(const double dt, const double fD)


\author Guido Montorsi
*/
class Fading_Channel  
{
public:

	//! Constructor and initialization of main parameters
	Fading_Channel(
		const double dt,			//!< Sampling interval
		const double fd=1.f,		//!< Doppler frequency
		const int maxdelay=1000		//!< Maximum delay of taps
		);

	virtual ~Fading_Channel();

	//! Add a tap to the channel
	void AddFadingTap(
		const double delaysin,//!< Delay of tap (number of samples, rounded to the closest integer)
		const double powersin //!< Power of the associated process in dB (or magnitude of coeff.)
		);

	//! Normalize the total power of processes 
	void Normalize(
		const double=1.	//!< Desired total power of processes (linear, default to one)
		);


	//! Set the Doppler frequency
	void SetDopplerFrequency(
		const double fD	//!< Desired Doppler frequency
		);

	//! Force the filter to a new random state.
	/*! All the content of the filters generating the taps coefficients are loaded
with new random values. The effect of a call to this method is that
of skipping a very large time interval or equivalently to generate a
new instance of the process. When the the Doppler frequency is zero,
a call to this method will generate a new set of constant random tap
coefficients.
	*/
	void Reset();

	//! Set the seed for the generation of all random processes 
	void SetSeed(const int seed);

	//! Simulation method
	/*! The output buffer can coincide with the input buffer */
	void Run(const int nsamples,		//!< Number of processed complex samples
		const double* input,    //!< Input complex samples
		double * output,		//!< Output complex samples
		double * fadingtaps=0	//!< Optional buffer to store fading taps		
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int nsamples,		//!< Number of processed complex samples
		const cmplx* input,    //!< Input complex samples
		cmplx* output,		//!< Output complex samples
		cmplx* fadingtaps=0	//!< Optional buffer to store fading taps		
		)
	{
	 Run(nsamples,		//!< Number of processed samples
		(const double*) input,    //!< Input complex samples
		(double*) output,		//!< Output complex samples
		(double*) fadingtaps	//!< Optional buffer to store fading taps		
		);
	}
#endif

	void Print(FILE* file = stdout);
	void ChangeFadingTaps(const double* delaysin,const double* powersin);
	void SetBirthDeath(const int rangein, const int periodin, const int step);
	void SetBirthDeath(const int range, const int period);
	void SetMoving(const double* Ain, const double *ffin);

	void Run_Frequency_Domain(
		const int N, // Number of frequencies
		const int nsamples,	// Number of blocks
		const double* input, //input  
		double * output,//output
		double* a=NULL);	// Channel transfer function

	double GetStrongestDelay(void);
	double MeasurePower();
	int Maxdelay();

	int* delays;			// Delays (number of chips)
	bool movingdelay;		// Moving delay on
	bool birthdeath;		// birthdeath on
	bool measure;

private:
	//! Set the Doppler frequency to zero. 
	void SetStatic(bool isstatic=true);
	int seed;				// Seed for the fading processes
	int	ntaps;				// number of taps of the fading filters
	double* h;				// Filter coefficients
	double **filters;		// Tapped delay lines for fading processes
	double* pointold;		// last fading values
	double* pointnew;			// future fading values


	double dt;				// Sampling interval of the input proces
	double fd;				// Dopppler frequency.

	/* Variable for moving */
	bool moving;
	double *A;
	double *B;
	double *ff;
	/* Variable for birth-death */
	int   range;			
	int   period;
	int	  step;
	int   pathdead;
	int	  maxdelay;			// Maxinum delay;

	int npaths;				// Number of paths
	double* levels;			// Levels of each path
	double *tapped;			// Tapped delay line for input process
	int oversampling;		// Oversampling of the fading process
	int time;				// time
	int pointer1;			// pointer for the fading processes TDL
	int pointer2;			// pointer for the input process TDL
	void NewSamples();

	double sigma;
	__int64	 nmeas;

	double* sen;
	bool isstatic;
	const static int alpha;



friend
Fading_Channel* Multipath(
		const double dtin, 
		const double fDin, 
		const int  npathsin, 
		const double* delaysin, 
		const double* powersin, 
		const int seedin);

 
friend Fading_Channel* Multipath_UMTSCase0(const double dt, const double fd);
friend Fading_Channel* Multipath_UMTSCase1(const double dt, const double fD);
friend Fading_Channel* Multipath_UMTSCase2(const double dt, const double fD);
friend Fading_Channel* Multipath_UMTSCase3(const double dt, const double fD);
friend Fading_Channel* Multipath_UMTSCase4(const double dt, const double fD);
friend Fading_Channel* Multipath_UMTSCase5(const double dt, const double fD);
friend Fading_Channel* Multipath_UMTSCase6(const double dt, const double fD);
friend Fading_Channel* Multipath_COST207_TU6( const double dt,const double fD);

};


Fading_Channel* Multipath(const double dtin,const double fDin,const int  npathsin,const double* delaysin,const double* powersin,const int seedin);
Fading_Channel* Multipath_UMTSCase0(const double dt, const double fd);
Fading_Channel* Multipath_UMTSCase1(const double dt, const double fD);
Fading_Channel* Multipath_UMTSCase2(const double dt, const double fD);
Fading_Channel* Multipath_UMTSCase3(const double dt, const double fD);
Fading_Channel* Multipath_UMTSCase4(const double dt, const double fD);
Fading_Channel* Multipath_UMTSCase5(const double dt, const double fD);
Fading_Channel* Multipath_UMTSCase6(const double dt, const double fD);

/*!\ingroup Interface
\brief Return the multipath fading channel TU6 specified in the project COST207
*/
Fading_Channel* Multipath_COST207_TU6(
									  const double dt, //!< Sampling interval
									  const double fD	//!> Doppler ferquency
									  );