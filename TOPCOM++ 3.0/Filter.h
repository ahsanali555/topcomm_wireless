#pragma once
#include <stdlib.h>
#include <stdio.h>
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class Filter and all its interfaces
*/
/*! Filter class and all its interface functions */
/*! \ingroup DSP
\brief Generic complex or real, FIR or IIR digital filter.


The user specifies:
- The number of feed-forward taps of the filter \f$N_z\f$
- The values of the feed-forward coefficients \f$z\f$
- The number of feedback taps (default to 0) \f$N_p\f$
- The values of the feed-feedback coefficients \f$p\f$
- A flag indicating if the coefficients are complex (default to NO)
- A flag indicating if the inputs are complex (default to NO).

Note that the first feedback tap is always neglected so that minimum number of 
active feedback taps to consider the
filter an IIR filter is 2. Also the correspondent value \f$p[0]\f$ is neglected.


The Run method has a behavior that depends on the above flags,
in particular the input vector is assumed to be complex if the
correspondent flag is on and the output is assumed to be complex if
the input is complex or the coefficients are complex. If the
pointer to the input vector is set to ``0" the filter inputs are
assumed to be all zero.

The method RunImpulse is used to run the filter over a sub-sampled sequence. 
It requires an additional parameter
\f$n_s\f$, such that \f$n_s-1\f$ represents the number of zeroes inserted between each 
sample of the input sequence. This
method is used for the shaping filter of the transmitter.

A number of interface functions for the filter are available to
construct filters for typical applications.
\see
- Concatenate()
- SRRC()
- SRJakes()
- Raised_Cosine()
- Exponential()
- GaussianFilter()
- Square()
- Echoes()
- OnePole()
- Filter_From_Mask()
- ComputeFIR()
- Chebyshev()
- Elliptic()
- Bessel()
- Butterworth()
- Truncate()
- Integrator()
- Derivator()
- IntDer()

\author Guido Montorsi
*/
class Filter
{
public:
	Filter();							
	virtual ~Filter();

	//! Set numerator and denominator of filter.
	void SetParameters(
		const int Nzin,				//!< Number of FF taps
		const double* forwin,		//!< Feed-forward taps
		const int Npin=0,			//!< Number of FB taps
		const double* backin=0,		//!< FB taps
		const bool hcmpx=false,		//!< Flag for complex coefficients.
		const bool icmpx=false		//!< Flag for complex inputs.
		);


	//! Reset Filter contents
	void Reset(
		const double* init = 0,		//!< Initial values of taps content
		const int period = 0		//!< Initial values of taps content
		);


	//! Run the Filter 
	/*! The output buffer can coincide with the input buffer. Overloaded with complex inputs. */
	void Run(
		int tics,				//!< Number of processed samples
		const double* input,	//!< Input samples	
		double* output			//!< Output samples
		);

	//! Generate a normalized (0 dB at DC) filter from its poles and zeroes
	void Set_Zeroes_and_Poles(
		const int Nzin,			//!< Number of zeroes
		const double* zeroes,	//!< Zeroes (Complex)
		const int Npin=0,		//!< Number of poles
		const double* poles=0,	//!< Poles (Complex)
		const bool hcmpx=false,	//!< Flag for complex  tap coefficients
		const bool icmpx=false	//!< Flag for complex inputs
		);	

	//! Scale the gain of filter.
	void SetGain(const double g		//!< Factor for gain.
		){gain*=g;};		

	//! Set Matched Filter (invert time axis for FIR filters)
	void SetMatched();

	//! Run the Filter over a impulse train
	void RunImpulse(int tics, //!< Number of processed samples.
		const double* input,  //!< Input samples.
		double* output,		  //!< Output samples 
		int ns				  //!< Number of inserted zero samples plus one.
		);

	//! Return the impulse response of the filter
	void Impulse(
		const int ntics,		//!< Number of samples to be generated.
		double* output			//!< Generated impulse response.
		);

	//! Return amplitude and phase of the transfer function of the filter	
	double* Amplitude_Phase(
		int m=14		//!< Base 2 logarithm of the number of considered samples of impulse response.
		);

	//! Return the energy of impulse response
	double Energy(
		const int l=1000	//!< Length of the considered impulse response.
		)
	{
		int i;
		double *temp= new double[l];
		Impulse(l,temp);
		double en=0.;
		for(i=0;i<l;i++)en+=temp[i]*temp[i];
		delete[] temp;
		return en;
	}
	//! Normalize the impulse response to have unitary energy.
	void Set_Unitary_Energy();

	//! Normalize the impulse response to have unitary DC.
	void Set_Unitary_DC(); 

	//! Print the Imuplse response to the desired output
	void PrintImpulse(
		const int tics,	//!< Number of generated samples of impulse response.
		FILE* file=stdout		//!< Output stream.
		);

	// Friends..
	friend class Gaussian_Process;

	double* line;		//!< line
	int Ntaps;			//!< Number of taps delay line
	double* forw;	    //!< Forward tap coefficients
	bool icmpx;			//!< Flag for complex inputs
	bool hcmpx;			//!< flag for complex coefficients

	//! Operator =			
	Filter& operator=(const Filter&);


	int Nz;				//!< Number of zeroes
	int Np;				//!< Number of poles
	int now;			//!< Pointer to delay line

#ifdef CTOPCOM
	//!
	void Run(
		int tics,				//!< Number of processed samples
		const cmplx* input,	//!< Input samples	
		cmplx* output			//!< Output samples
		)
	{
		Run(
			tics,				
			(const double*) input,	
			(double*) output			
			);
	}

	void RunImpulse(int tics, //!< Number of processed samples.
		const cmplx* input,  //!< Input samples.
		cmplx* output,		  //!< Output samples 
		int ns				  //!< Number of inserted zero samples plus one.
		)
	{
		RunImpulse(tics, 
			(const double*) input, 
			(double*) output,		  
			ns				  
			);
	}
#endif


private:
	void RunCC(int tics,const double* input, double* output);
	void RunCR(int tics,const double* input, double* output);
	void RunRC(int tics,const double* input, double* output);
	void RunRR(int tics,const double* input, double* output);

	bool matched;	//!< Matched filter (inverted time)



	bool    IIR;		//!< Flag for IIR type filter
	double* back;		//!< Backward tap coeficients
	double gain;

	   /* Interface functions */
    friend
    Filter* Filter_From_Mask(const int N, 
            const double *f,
            const double *H,
            const double *F,
            const double ds, 
            const double perc, 			  
            const int interp,
            const bool hcmplx, 
            const bool icmplx 
            );


    friend
    Filter* ComputeFIR(const int N1, 
            const int N, 
            const double *H, 
            const double *F, 
            const double perc,
            const int off,
            const int window, 
            const bool hcmplx,
            const bool icmplx 
            );



    friend
    Filter* Concatenate(Filter* a, Filter* b, const int N);



    friend
    Filter* Truncate(
            Filter* ref, 
            const double bandwidth, 
            const double rolloff
            );


    friend
    Filter* SRRC(const double alpha, 
            const int ns, 
            const int N, 
            const bool iscmplx, 
            const double delfrac 
            );

    friend
    Filter* DVBRCS(const double alpha, 
            const int ns, 
            const int N,
            const bool iscmplx, 
            const double BWnorm
            );


    friend
    Filter* SRJakes(
            const int over, 
            const bool iscmplx,
            const double eps
            );


    friend
    Filter* Rounded(
            const int over, 
            const bool iscmplx,
            const double eps 
            );


    friend
    Filter* Raised_Cosine(const double alpha,
            const double gamma, 
            const double dt,
            const int N, 
            const bool iscmpx
            );


    friend
    Filter* Time_RC(const double alpha, 
            const double gamma, 
            const int ns,
            const bool iscmpx
            );

 
    friend
    Filter* Exponential(const double decay,
            int N, 
            const bool iscmplx
            );


    friend
    Filter* GaussianFilter(const double sigma, 
            int N, 
            const bool iscmplx
            );



    friend
    Filter* Square(const double W,
            int N, 
            const bool iscmplx
            );

    friend
    Filter* OnePole(const double alpha, bool icmplx);

    friend
    Filter* Butterworth(
            const double fc,
            const int N, 
            const bool iscmplx
            );


    friend
    Filter* Chebyshev(
            const double ripdb, 
            const double fc,
            const int order, 
            const bool iscmplx
            );


    friend
    Filter* Elliptic(
            const double ripdb,
            double fc, 
            double fa, 
            const int order, 
            const bool iscmplx
            );



    friend
    Filter* Bessel(
            const double fc, 
            const int n, 
            const bool iscmplx
            );


    friend
    Filter* IntDer(int der,  
            const int degree,
            bool iscmplx 
            );

  
    friend
    Filter* Derivator(
            const int degree , 
            bool iscmplx	 
            );


    friend
    Filter* Integrator(
            const int degree , 
            bool iscmplx  
            );

    friend
    Filter* FIR_Estimation_Filter(
            Filter* Sd,
            Filter* Se, 
            const int length, 
            const int del, 
            const bool icmplx
            );

    friend
    
    Filter* Hilbert_Raised_Cosine(const double alpha, const double gamma, const double dt, const int N1);

};

/* Interface functions */

/*! \ingroup Interface
\brief Return a filter with a given mask, specified through a set of points to be interpolated
 */

Filter* Filter_From_Mask(const int N, //!< Number of points
        const double *f, //!< Frequencies
        const double *H, //!< Amplitudes [dB]
        const double *F, //!< Phases [degree]
        const double ds, //!< Spacing of interpolate samples
        const double perc = 1., //!< Fraction of kept energy in the filter				  
        const int interp = 0, //!< Interpolation method.
        const bool hcmplx = false, //!< The coefficient of filter are complex?
        const bool icmplx = false //!< The inputs of filter are complex ?
        );

/*! \ingroup Interface
\brief Return a filter with a given mask, specified through a set of equally spaced points
 */

Filter* ComputeFIR(const int N1, //!< Number of frequency samples
        const int N, //!< Number of provided frequency samples
        const double *H, //!< Amplitude of samples 
        const double *F, //!< Phase of samples 
        const double perc = 1., //!< Percentage of energy for truncation
        const int off = 0, //!< Offset of first sample.
        const int window = 1, //!< Type of truncating window.
        const bool hcmplx = false, //!< The coefficient are complex?
        const bool icmplx = false //!< The inputs are complex ?
        );


/*! \ingroup Interface
\brief Return a filter obtained serially concatenating two filters. 
 */

Filter* Concatenate(Filter* a, Filter* b, const int N);


/*!\ingroup Interface 
\brief Return a Filter with transfer function truncated into a given bandwidth 

 */

Filter* Truncate(
        Filter* ref, //!< Starting filter
        const double bandwidth, //!< Bandwidth (smaller than one)
        const double rolloff = 0 //!< Roll-off of truncating window
        );

/*! \ingroup Interface
\brief Returns a  square root raised cosine filter.

Close formula for the impulse response is used. 
The length of the filter, if not specified by the user, is computed 
automatically by suitably truncating the impulse response where it becomes 
negligible.

\see Raised_Cosine()
 */

Filter* SRRC(const double alpha, //!< Roll-off.
        const int ns, //!< Number of samples per symbol.
        const int N = -1, //!< Length of generated FIR (default to automatic length computation).
        const bool iscmplx = false, //!< Flag for complex inputs.
        const double delfrac = 0. //!< fractional delay .
        );


Filter* DVBRCS(const double alpha, //!< Roll-off.
        const int ns, //!< Number of samples per symbol.
        const int N = -1, //!< Length of generated FIR (default to automatic length computation).
        const bool iscmplx = false, //!< Flag for complex inputs.
        const double BWnorm = 1 //!< Normalized bandwidth
        );

/*! \ingroup Interface
\brief Returns a  square root Jake's shaped FIR filter.

The FIR filter is obtained by inverse DFT and truncation of the function
\f[ 
H(f)=\sqrt{\frac{1}{\sqrt{1-(f/f_D)^2}}}.
\f]
 */

Filter* SRJakes(
        const int over, //!< Oversampling factor fd==1
        const bool iscmplx, //!< Flag for complex inputs.
        const double eps = 0.01 //!< Error for truncation
        );

/*! \ingroup Interface
\brief Return a  Rounded filter (for Doppler Spectra) 
 */

Filter* Rounded(
        const int over, //!< Oversampling factor fd==1
        const bool iscmplx, //!< Flag for complex inputs.
        const double eps = 0.01 //!< Error for truncation
        );

/*! \ingroup Interface
\brief Return a  raised cosine filter raised to an arbitrary power.

The filter is computed by inverse DFT of the frequency response.
\see SRRC()
 */

Filter* Raised_Cosine(const double alpha, //!< Roll-off
        const double gamma, //!< Exponent of the Raised Cosine
        const double dt, //!< Sampling interval
        const int N, //!< Length of generated FIR (must be odd)
        const bool iscmpx = false //!< Flag for complex inputs
        );

/*! \ingroup Interface
\brief Return a  filter with raised cosine IMPULSE RESPONSE.

See Raised_Cosine() or SRRC() for more common filter with raised cosine transfer function.

\see SRRC()
\see Raised_Cosine()
 */

Filter* Time_RC(const double alpha, //!< Roll-off
        const double gamma, //!< Exponent of the Raised Cosine
        const int ns, //!< Samples per symbol
        const bool iscmpx = false //!< Flag for complex inputs
        );

/*! \ingroup Interface
\brief Return a  FIR filter with (truncated) exponential impulse response
 */

Filter* Exponential(const double decay, //!< Time constant.
        int N = 0, //!< Length of generated FIR (0=Automatic computation).
        const bool iscmplx = false //!< Flag for complex inputs.
        );

/*! \ingroup Interface
\brief Return a  FIR filter with (truncated) Gaussian impulse response
 */

Filter* GaussianFilter(const double sigma, //!< Variance.
        int N = 0, //!< Length of generated FIR (0=Automatic computation).
        const bool iscmplx = false //!< Flag for complex inputs.
        );


/*! \ingroup Interface
\brief Return a  FIR filter with square transfer function.
 */

Filter* Square(const double W, //!< Bandwidth.
        int N = 0, //!< Length of generated FIR (0=Automatic computation).
        const bool iscmplx = false //!< Flag for complex inputs.
        );
/*! \ingroup Interface

  \brief Return an IIR filter with impulse response \f$h_n=\alpha^n\;\;\;n=0,\dots,\infty\f$
 */

Filter* OnePole(const double alpha, bool icmplx = false);

/*!\ingroup Interface
\brief Return a Butterworth filter.

User specifies the filter features by setting the cut-off frequency and the filter order \f$N\f$.
The interface uses the interface Chebyshev()
 */

Filter* Butterworth(
        const double fc, //!< Cut-off frequency.
        const int N, //!< Filter order.
        const bool iscmplx = false //!< The inputs are complex ?
        );

/*!\ingroup Interface
\brief Return a Chebyshev filter.

User specifies the filter features by setting the ripple, the cut-off
frequency and the filter order \f$N\f$.
 */

Filter* Chebyshev(
        const double ripdb, //!< Ripple in dB.
        const double fc, //!< Cut-off frequency.
        const int order, //!< Filter order.
        const bool iscmplx = false //!< The inputs are complex ?
        );

/*!\ingroup Interface
\brief Return an Elliptic filter.

User specifies the filter features by setting the ripple in passband, the passband frequency, 
the stopband frequency and the filter order \f$N\f$.
 */

Filter* Elliptic(
        const double ripdb, //!< Ripple in dB.
        double fc, //!< Passband frequency.
        double fa, //!< Stopband frequency.
        const int order, //!< Filter order.
        const bool iscmplx = false //!< The inputs are complex ?
        );



/*!\ingroup Interface
\brief Return a Bessel filter.

User specifies the filter features by setting the cutoff frequency, 
and the filter order \f$N\f$.

 */

Filter* Bessel(
        const double fc, //!< Cut off frequency. 
        const int n, //!< Filter order.
        const bool iscmplx = false //!< The inputs are complex ?
        );

/*!\ingroup Interface
\brief Return a filter that performs integration or derivation of an input signal.

The operations are performed over an interpolated signal from the input sequence.
The user can specify the desired degree of the interpolating polynomial, which by default is 1, 
meaning that linear interpolation is performed.
The delay introduced by the block is half the degree of the interpolating polynomial
The user can specify the desired degree of the interpolating polynomial, which by default is 3.

 */

Filter* IntDer(int der, //!< >0 order of derivative <0 Order of integration 
        const int degree = 1, //!< Degree of the polynomial (default to cubic interpol)
        bool iscmplx = false //!< Is the signal complex?	 
        );

/*!\ingroup Interface
\brief Return a  filter that performs derivation of the input signal.

The derivation is performed over an interpolated signal from the input sequence.
The user can specify the desired degree of the interpolating polynomial, which by default is 1, 
meaning that linear interpolation is performed.
The delay introduced by the block is half the degree of the interpolating polynomial

 */

Filter* Derivator(
        const int degree = 1, //!< Degree of the polynomial (default to linear interpol)
        bool iscmplx = false //!< Is the input signal complex?	 
        );

/*!\ingroup Interface
\brief Return a filter that performs integration of the input signal.

The integration is performed over an interpolated signal from the input sequence.
The user can specify the desired degree of the interpolating polynomial, which by default is 1, 
meaning that linear interpolation is performed.
The delay introduced by the block is half the degree of the interpolating polynomial

 */

Filter* Integrator(
        const int degree = 1, //!< Degree of the polynomial (default to linear interpol)
        bool iscmplx = false //!< Is the signal complex?	 
        );

/*! \ingroup Interface
\brief Return an optimal non causal Wiener filter for optimal estimation
of data in uncorrelated noise.

The filter has real tranfer function
\f[
H(z)=\frac{N_d(z)\cdot D_n(z)}{N_n(z)\cdot D_d(z)+N_d(z)\cdot D_n(z)}
\f]


where \f$N_d(z), D_n(z), N_n(z), D_d(z)\f$ 
are the squared absolute values of the numerator (N) and denominator
(D) of the data (d) and noise (n) filters 
 */

Filter* FIR_Estimation_Filter(
        Filter* Sd, //!< Data filter
        Filter* Se, //!< Noise filter.
        const int length, //!< Length of filter
        const int del = 0, //!< Default to filtering
        const bool icmplx = false //!< The input  are complex ?
        );


Filter* Hilbert_Raised_Cosine(const double alpha, const double gamma, const double dt, const int N1 = 200);
