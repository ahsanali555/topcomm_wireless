// SCTCM_Encoder.h: interface for the SCTCM class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "TCM_Encoder.h"
#include "TCM_Decoder.h"
#include "Trellis_Encoder.h"
#include "Interleaver.h"
#include "SISO_Decoder.h"
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of class SCTCM_Encoder 
*/
/*! \ingroup Modems 
\brief Serially concatenated trellis coded modulation encoder.

The encoder is composed by an outer binary convolutional encoder (Trellis_Encoder),
a binary block interleaver (Interleaver) and an inner trellis coded modulation encoder (TCM_Encoder). 
The SetParameters() method of this class requires pointers to the three constituent blocks, 
that have to be defined and initialized in the main program.


The SCTCM encoder operates in terminated mode, i.e., both outer and inner encoder are terminated 
and the interleaver operates in block mode so that the state of the SCTCM encoder at the end of each encoded block 
is fixed and the encoder as a whole is a block encoder.

The variable \f$K\f$ and \f$N\f$ are public so that the user can access them to properly 
initialize the size of input and
output buffers. The interleaver performs the permutation at bit level.

A single tic of the Run() method corresponds to the encoding of a whole coded block and thus requires 
\f$K\f$ input bits
and provide a set of \f$N\f$ complex or real values representing the constellation points.

For an example of its use see e.g. the test program "test_SCTCM.cpp".

\author Guido Montorsi

*/
class SCTCM_Encoder  
{
public:
	SCTCM_Encoder();
	virtual ~SCTCM_Encoder();

	//! Set the main parameters of encoder.
	void SetParameters(Trellis_Encoder* outer,	//!< Outer trellis encoder.
		const Interleaver*  interleaver,		//!< Binary interleaver.
		TCM_Encoder* inner						//!< Inner TCM encoder.
		);				

	//! Simulation method
	void Run(const int tics,					//!< Number of encoded blocks.
		const int* inp,							//!< Input bits.
		double* out								//!< Output constellation points.
		);	

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,	//!< Number of encoded blocks
		 const int* inp,		//!< Input bits
		 cmplx* out				//!< Complex Output contellation points
		 )
	{
	Run(tics,inp,(double*) out);
	return;
	}
#endif

	friend class SCTCM_Decoder;
	int K;			//!< Number of input bits.
	int I;			//!< Interleaver length.
	int N;			//!< Number of output consetllation points.

private:
	Trellis_Encoder *outer;
	TCM_Encoder   *inner;
	const int* perm;
	int* buff;
};

/*! \file */
/*! \ingroup Modems 
\brief Serially concatenated trellis coded modulation decoder.

The SCTCM_Decoder is the block associated to the SCTCM_Encoder that performs iterative decoding for serially
concatenated trellis coded modulations. 

The parameters for the construction of SCTCM decoder are the pointer to the corresponding encoder,
that determines the encoding scheme structure, and the number of
iterations required for the decoding. Optional parameters are the
window for the initialization and updating of the SISO algorithm
\f$T_W\f$ and \f$U_W\f$, and  the precision factor \f$f\f$ used for the
representation of the input metrics.

The input to the Run() method are complex or real samples of the received signal. 
The samples are assumed to be
normalized with respect to a constellation with unitary energy.

A single tic of the Run() method corresponds to the decoding of a whole coded block 
and thus requires \f$N\f$ input
samples and provide a set of \f$K\f$ decoded information bits. The
Run method takes as input the \f$I,Q\f$ samples of the received signals and fed them to a Demodulator
 that computes their LLRs with the RunSoft() method.

The following block is an iterative SCCC decoder  
that uses as inner code the trellis embedded in the TCM scheme.

For an example of its use see e.g. the test program "test_SCTCM.cpp".

\author Guido Montorsi
*/
class SCTCM_Decoder
{
public:
	SCTCM_Decoder();
	virtual ~SCTCM_Decoder();

	//! Set the main parameters of encoder.
	void SetParameters(
		const SCTCM_Encoder* e,		//!< Pointer to the correspondent encoder
		const int niter,			//!< Number of iterations.
		const int nini     =50,		//!< Training window of SISO processors.
		const int ngroup   =100,	//!< Updating window of SISO processors.
		const double factor=8.		//!< Scaling factor of input metrics.
					  );
	//! Set the value of the variance of noise and the precision of the llr representation
	void SetSigma(const double sigma2,		//!< Variance of noise 
				const double factin=8.		//!< Conversion factor \f$f\f$ between actual llrs values \f$\lambda_i\f$ (double) and their quantized version (integer).
				);

	//! Change the number of iterations of decoder.
	void TuneParameters(const int niterin	//!< Number of iterations.
		);

	//! Run the decoder
	void Run(
			const int tics,		//!< Number of decoded blocks
			const double* inp,	//!< Input received complex samples.
			int* out			//!< Output (soft) bits.
	);	

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics,	//!< Number of decoded blocks
		 const cmplx* inp,		//!< Input received complex samples
		 int* out				//!< output (soft) bits
		 )
	{
	Run(tics,(double*) inp,out);
	return;
	}
#endif

	//! Turn on and off the stopping rule.
	void SetStop(const bool stopin=true){stop=stopin;}

	Demodulator *dem;		//!< Embedded Demodulator. 
	SISO_Decoder * iSISO;	//!< Embedded Inner SISO.
	SISO_Decoder * oSISO;	//!< Embedded Outer SISO.

private:
	const SCTCM_Encoder* refcod;
	int niter;
	int *llr;		// Temporary llr
	int *ext1;		// Temporary llr
	int *ext2;		// Temporary llr
	int nti;		// Number of trellis steps inner
	int nto;		// Number of trellis steps outer
	int nsymbols;
	bool stop;
};
