// Demodulator.h: interface for the Demodulator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Modulator.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of class Demodulator*/
/*! \ingroup Modems
\brief Generic one dimensional or two dimensional demodulator.

User can specify the constellation and the mapping between binary input and
one or two dimensional constellation points.
The module provides both hard decoding and soft decoding.

Hard demodulation: the demodulator implements a minimum-distance decision rule, 
corresponding to optimum maximum-likelihood demodulation.
Depending on the value of a flag set by the user, 
the output vector may contain either the symbols or the binary digits.

Soft demodulation: the demodualtor generates a vector containing 
the log-likelihood ratios \f$\lambda_i\f$ of the \f$M\f$ possible transmitted waveforms, 

For an example of its use see e.g. the test program "test_PCCC.cpp".

\author Guido Montorsi
*/
class Demodulator  
{
public:
	Demodulator();
	virtual ~Demodulator();
	//! Set the main parameters of the demodulator
	void SetParameters(const Modulator* mod //!< Modulator
						);
	//! Set the value of the variance of noise and the precision of the llr representation
	void SetSigma(const double sigma2,			//!< Variance of noise 
				const double factin=8.		//!< Conversion factor \f$f\f$ between actual llrs values \f$\lambda_i\f$ (double) and their quantized version (integer).
				);

	//! Set quantization of innput values
	void SetQuant(const int nbit, const double scale);

	//! Hard decoding demodulator 
	void Run(
		const int ntics,				//!< Number of input symbols
		const double* inp,				//!< Input signal 
		int* out,						//!< Output symbols or bits
		const bool bin=true				//!< If bin is true, the bits are output instead of symbols
		);

	//! Soft decoding demodulator 
	void RunSoft(
		const int tics,					//!< Number of input symbols
		const double* inp,				//!< Input signal 
		int* out,						//!< Output LLRs on symbols or bits 
		const bool bin=true	,			//!< If bin is true, the LLRs of bits are output instead of those of symbols
		const double* chabs=0			//!< Scaling factor for the variance of noise (default to 1)
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	//! Hard decoding demodulator (complex input)
	void Run(
		const int ntics,				//!< Number of input symbols
		const cmplx* inp,				//!< Complex Input signal 
		int* out,						//!< Output symbols or bits
		const bool bin=true				//!< If bin is true, the bits are output instead of symbols
		)
	{
	Run(ntics,(const double*)inp,out, bin);
	return;
	}

	//! Soft decoding demodulator (complex input)
	void RunSoft(
		const int tics,					//!< Number of input symbols
		const cmplx* inp,				//!< Complex Input signal 
		int* out,						//!< Output LLRs on symbols or bits 
		const bool bin=true	,			//!< If bin is true, the LLRs of bits are output instead of those of symbols
		const cmplx* chabs=0			//!< Scaling factor for the variance of noise (default to 1)
		)
	{
		RunSoft(tics,(const double*) inp,out,bin,(const double*)chabs);
		return;
	}

#endif

	/** \cond EEE */
	//! Soft decoding demodulator for observation affected by phase jitter
	/*! Computes exact quantized LLR of symbols when the received samples are affected by a
	random Tichonov phase shift. Requires that the phase jitter variance is passed to the block
	through the method ::SetSigmaTheta. Use RunSoft otherwise.
	*/
	void RunSoftTheta(
		const int tics,					//!< Number of input symbols
		const double* inp,				//!< Input signal 
		int* out,						//!< Output LLRs on symbols or bits 
		const bool bin=true			//!< If bin is true, the LLRs of bits are output instead of those of symbols
		);

	//! Run pragmatic demodulator with out LLR
	void RunGraySubopt(const int tics, 
				  const double* inp, 
				  int* out
						  );

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void RunSoftTheta(
		const int tics,					//!< Number of input symbols
		const cmplx* inp,				//!< complex Input signal 
		int* out,						//!< Output LLRs on symbols or bits 
		const bool bin=true			//!< If bin is true, the LLRs of bits are output instead of those of symbols
		)
	{
		RunSoftTheta(tics,(const double*) inp,out,bin);
		return;
	}

	//! Run pragmatic demodulator with out LLR
	void RunGraySubopt(const int tics, 
				  const cmplx* inp, 
				  int* out
	)
	{
	RunGraySubopt(tics,(const double*) inp, out);
	return;
	}
#endif

	//! Set the value of the ratio between the variances of phase jitter and additive noise
	void SetSigmaTheta(
		const double K	//!< Ratio between variances of phase and noise 
		);

	//! Compute an estimate of the transmitted signal from apriori on the bits or symbols
	void Estimate(
			const int tics,				//!< [in] Number of symbols
			const int* inp,				//!< [in] Input LLRs on symbols or bits 
			double* out,				//!< [out] Output estimated signal 
			double* outvar,				//!< [out] Variance of estimate
			const bool bin=true			//!< [in] If bin is true, the LLRs of bits are provided instead of those of symbols
			);
	/** \endcond */


	int  maxtab;					//!< Size of max* lookuptable
	int* ic;						//!< max* lookup table

	friend class Tuning_Centroids;

private:
	//! Initialize output of soft binary LLR
	void SetSoftBin(); 
	bool soft;						//!< Is soft decoding?
	const Modulator * refmod;		//!< Modulator
	double fact;					//!< Conversion factor
	double a;						//!< a=fact/2/sigma2
	int* LF;						//!< Likelihood functions
	int* split;
	int	  ni;						//!< Number of output alphabets
	int   nli;						//!< Size of output vector
	double sigma2;					//!< noise variance (for exact LLR computation)
	double K;						//!< Ratio between noise and phase noise variances
	double* bestt;					//!< Best rotations for demodulator with phase noise
	friend class DVBSX_Receiver;
};
