// Spectrum.h: interface for the Spectrum class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdio.h>
#include "FFT.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file  
\brief Declaration of the class Spectrum
*/
/*! \ingroup Measurements 
\brief  Evaluate the power spectrum of a discrete process.

 The Spectrum class implements a power spectral density (PSD) probe, using the
 Bartlett periodogram method, based onthe following algorithm:
- Divide the entire input sequence \f$x[n]\f$ into \f$K\f$ subsequences \f$x_k[n]\f$ 
(\f$k=0,...,K-1\f$) of length \f$N=2^n\f$ samples each.
- For each subsequence, evaluate the \e Periodogram \e \f$\hat{P}^{K}_x[i]\f$ using one of the 
following windowing signals: Hamming window, Hanning window or non window.

The user specifies:
- The \f$\log_2\f$ of the number of samples to perform the FFT.
- The sampling interval of the input signal.
- The name of output file in which the histogram is stored.
- The type of windowing (Hamming, Hanning or none).

\author Guido Montorsi
*/

class Spectrum  
{
public:
	Spectrum();
	//! Constructor with parameters
    Spectrum(int nsampin,	//!< log2(# of samples) to perform the fft
		double T,			//!< Sampling interval
		char* name,			//!< Name of output file
		int windtype=1	//!< Window of data: 0=non window, 1=Hamming, 2=Hanning (default: Hamming)
	);
	
	virtual ~Spectrum();

	//! Set the main parametrs of the block
	void SetParameters(int nsampin,	//!< log2 Samples to perform the fft
		double T,			//!< Sampling interval
		const char* name,			//!< Name of output file
		int windtype=1	//!< Window of data: 0=non window, 1=Hamming, 2=Hanning (default: Hamming)
		);

	//! Reset statistic
	void Reset();

	//! Run the block 
	/** First overload: real and imaginary parts are in two separate vectors */
	void Run(const int ntics,//!< Number of input samples 
		const double * Real,		//!< Real part of input signal 
		const double * Imag	//!< Imaginary part of input signal 
		);
	//! Run the block 
	/** Second overload: real and imaginary parts are in a single vector in even and odd positions */
	void Run(const int tics,	//!< Number of input samples 
			const double* Complex	//!< Input signal (real and imaginary parts are in a single vector in even and odd positions)
			);	
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,   //!< Number of input samples
			 const cmplx* Complex		//!< Compelx Input signal 
			 )
	   {
	   Run(tics,(double*) Complex);
	   }
#endif


	//! Print the spectrum in the output file
	void Print();

	//! This flag allows to print the autocorrelation instead of the spectrum
	bool autocorr;

	//! Evaluate the bandwidth that contains a given percentage of the total power
	double Bandwidth(double percentage //!< desired percentage
		);
	double* buffer;	//!< Vector containing the periodogram (or autocorrelation)

private:
	double norm;	//!< Window normalization factor (the energy must be equal to 1)
	FFT a;	//!< FFT member
	int windtype;	//!< Type of windowing
	int nsamp;		//!< Number of samples
	int nn;			//!< log2 of numnber of samples to perform the fft
	int kp;			//!< Time index
	int nblock;		//!< Number of periodograms
	double* fact;	//!< Vector containing the window samples
	double* aver;	//!< Average of periodograms
	FILE* file;		//!< Name of output file
	double T;		//!< Sampling interval
};
