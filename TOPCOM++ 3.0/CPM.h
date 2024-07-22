/*!\file
\brief Declaration of CPM_Encoder and CPM_Decoder class and their interface functions.
*/

#pragma once

const double pi = 3.1415926535897932384626433832795;

enum CPM_type {REC, RC, GAUSSIAN, AV}; 

#include "Trellis.h"
#include "ShiftFrequency.h"
#include "Viterbi.h"
#include "Filter.h"
#include "Delay.h"
#include "EyeDiagram.h"
/*! \ingroup Modems
\brief Generic CPM modulator.
\author Guido Montorsi

 The user can set 
	- The pulse type among REC, raised cosine and Gaussian.
	- the modulation index in the form h=k/p 
	- The pulse length 
	- The mapping of bits to CPM levels

 The modulator implementation is based on the CPM representation described in [Bixio] and
 so the generate signal have a tilted frequency, to change it use 
 CPM_Modulator::Tilt_Frequency.
 Terminated mode is supported (see CPM_Modulator::Set_Terminated).
*/

class CPM_Modulator
{
public:
	CPM_Modulator();
	//! Constructor with parameters, see CPM_Modulator::SetParameters 	
	CPM_Modulator(int n,int K,int P, int L, CPM_type type,int nsamples,bool gray=false);
	
	//! Set the parameters of the CPM modulator
	void SetParameters(		
			int n,			//!< Number of bits per modulation signal.
			int K,			//!< Numerator of modulation index.
			int P,			//!< Denominator of modulation index.
			int L,			//!< Length of frequency pulse
			CPM_type type,	//!< Type of frequency pulse.
			int nsamples,	//!< Number of samples generated 
			bool gray=false	//!< Type of mapping of bits into CPM levels.
		);	
	//! Set the parametrs of a GMSK modualator
	/*! Method to be used in alternative to SetParameters */
	void SetGMSK(const double BT,	//!< Product \f$BT\f$ that usually defines the GMSK modualtion scheme
		const int nsamplesin		//!< Number of samples generated
		);

	//! Set the tilt of the frequency of generated signal
	void Tilt_Frequency(	
		const bool tilt=true//!< Flag to indicate if the tilt is desired or not.	
	);					

	virtual ~CPM_Modulator();
	//! Reset the state of the encoder	
	void Reset(){state=0;if(!tilt)shiftf->tempo=0;};

	//! Set the modulator to operate in terminated mode.
	/*! When this method is called the behavior of the CPM_Modulator::Run method changes.
	The state of the encoder is	periodically forced back to the zero state. 
	In this case one tics of the Run method generates
	a sequence of signals. The number of input bits required depends on the number
	of terminating steps required.*/
	void SetTerminated(
		int nsymbols,		//!< Number of trellis steps for termination
		bool before=true	//!< Before or after termination?
		);

	void Run(int tics,		//!< Number of tics 
		const int* bits,	//!< Input bits
		double* output,		//!< Output signal
		int * sigindex=0		//!< Output signal index
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(int tics, //!< Number of symbols 
		const int* bits,//!< Input bits
		cmplx* output, //!< complex Input samples
		int *sigindex=0			//!< Output LLRs.
		)
	{
		Run(tics,bits,(double*)output,sigindex);
		return;
	}
#endif
	Trellis* trellis;	//!< Trellis of encoder.
	int nsamples;		//!< Number of samples per signal.
	int terminated;		//!< Flag that indicate if the encoder is terminated.
	int nsymbols;		//!< For terminated CPM number of CPM symbols in block.
	int nterm;			//!< Required number of  steps to force the encoder to zero.
	CPM_type	type;   //!< Type of frequency pulse
	int K;				//!< Numerator of modulation index.
	int P;				//!< Denominator of modulation index.
	int L;				//!< Length of frequency pulse.
	int M;				//!< Cardinality of CPM signal.
	double h;			//!< Modulation index.
	int ncpm;			//!< Number of bits.
	double alpha_RC;	//!< Phase response of the weighted average (AV) CPM pulse shape used in DVBRCS2
	

	void SetFeedback(const int type, const int *feedback);

	friend class CPM_Demodulator;
	friend class SCCPM_Encoder;
	friend class SCCPM_Decoder;

private:
	void COMPUTE_TRELLIS_AND_SIGNALS();

	int* termination;	// 
	void Signal(int symbol, double* Real, double* Imag, int tau) const;
/* parameters */
	double (*q)(double,int);  //   
	bool tilt;
	ShiftFrequency * shiftf;	// For unbiased central frequency

	int u;      // input symbol
	int state;     // current state
	double **waver;  // Signals (Real)
	double **wavei;  // Signals (Imag)
//	double *Real;    // Output SIgnal
//	double *Imag;    // 
    bool gray;
	double *interfc;	// Phase pulse for the GMSK
};

//CPM_Modulator* GMSK(const float BT);


/*! \ingroup Modems
\brief Generic CPM memoryless demodulator  (see also CPM_Decoder) 
\author Guido Montorsi

This class is used to compute only the metrics associated to the memoryless modulator
of the CPM scheme. The class is used in iterative decoders (SCTCM) where the convolutional
structure  of CPM is used as inner code in a serial concatenated scheme.
Use class CPM_Decoder for a full CPM Decoder.
*/
class CPM_Demodulator  
{
public:
	CPM_Demodulator();
	//! Constructor with parameters, see CPM_Demdulator::SetParameters 	
	CPM_Demodulator(const CPM_Modulator*, const int nfilters);
	virtual ~CPM_Demodulator();
 	//! Set the parameters of the CPM demodulator
	void SetParameters(const CPM_Modulator*,		//!< Pointer to the CPM modualator 
		const int nfilters=0						//!< Number of matched filters
		);
   //! Run the demodulator
	void Run(const int tics,	//!< Number of symbols 
		const double* input,	//!< Input samples
		int *llr				//!< Output LLRs
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics, //!< Number of symbols 
		const cmplx* input, //!< complex Input samples
		int *llr			//!< Output LLRs.
		)
	{
		Run(tics,(const double*)input,llr);
		return;
	}
#endif

	//! Set the value of \f$E_s/N_0\f$ necessary to correclty compute the LLRs
	void SetNoise(const double esn0		//!< Signal to noise ratio
		);

	int Neff;// Total number of signals

	friend class CPM_Decoder;
	friend class SCCPM_Decoder;
	friend class DVBRCS_Receiver;

private:
	void ComputeFilters(const int nfilters);

/* Input parameters */
	int nfilters;		// Number of pseudo-matched filters
	int nsamples;		// Number of samples per symbol

/* derived parameters */
//	double (*q)(double,int);  // Phase pulse
	double h;			// Modulation index

/* inner memory*/
	double *r1r;		// Output of matched filters (Re)
	double *r1i;		// Output of matched filters (Im)
	double **hr;		// Impulse response of pseudo-MF (Real part) 
	double **hi;		// Impulse response of pseudo-MF (Imaginary part)
	double **mr;		// Signal vector in the pseudo-space (Re)
	double **mi;		// Signal vector int the pseudo-space (Im)
	int nsignals;		// Counter of received signals
	double *energy;      // Energy at the output of each matched filter
//	double *Lambdar;	// Covariance matrix (Re)
//	double *Lambdai;	// Covariance matrix (Im)
	double factor;		// N0 Varaince of noise 
	FILE* file;	
	const CPM_Modulator *mod; // Paired modulator
};


/*! \ingroup Modems
\brief Optimal (ML) CPM demodulator and decoder.
\author Guido Montorsi
  This decoder performs Maximum-Likelihood sequence decoding of CPM signals. The
  received signals is filtered through a minimal set of filters that provide the sufficient statistics 
  and then ML sequence estimation is performed through the Viterbi algorithm.

*/
class CPM_Decoder  
{
public:
	CPM_Decoder();
	virtual ~CPM_Decoder();
	//! Set the main parameters of the decoder
	void SetParameters(const CPM_Modulator* mod,	//!< CPM reference modulator.
					   const int nfilters,			//!< Maximum number of filters in the front-end
					   const int delay,				//!< Window or latency of the Viterby algorithm
					   const int ini=10			//!< Updating window of Viterbi algorithm.
);
	//! Run the decoder
	void Run(const int tics, //!< Number of symbols 
		const double* input, //!< Input samples
		int *llr			//!< Output LLRs.
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics, //!< Number of symbols 
		const cmplx* input, //!< complex Input samples
		int *llr			//!< Output LLRs.
		)
	{
		Run(tics,(const double*)input,llr);
		return;
	}
#endif

private:
	CPM_Demodulator *dem;
	Viterbi			*dec;
	int *llr;		// Temporary llr
	int* temp;
};

/*! \ingroup Modems
\brief Simple linear receiver for Gaussian MSK
\author Guido Montorsi
*/

class Simple_GMSK_Decoder
{
public:
	Simple_GMSK_Decoder();
	virtual ~Simple_GMSK_Decoder();

	//! Run the receiver
	void Run(
		int tics,				//!< Number of symbols 
		const double* input,	//!< Input samples (8 times the number of symbols)
		int *llr				//!< Output LLRs.
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int tics, //!< Number of symbols 
		const cmplx* input, //!< complex Input samples
		int *llr			//!< Output LLRs.
		)
	{
		Run(tics,(const double*)input,llr);
		return;
	}
#endif
	//! Reset the receiver
	void Reset(const int off=0);

private:
	Filter* Equal;
	Filter* C0;
	Delay * delay;
	double* buff;
	int off;
//	EyeDiagram *eye;
};

double qRC(double,int);
double qREC(double,int);
double qAV(double,int,double);

