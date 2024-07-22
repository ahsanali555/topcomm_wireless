// AWGN_Channel.h: interface for the AWGN_Channel class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdlib.h>
#include <math.h>

/*! \file 
\brief Declaration of the class AWGN_Channel 
*/

/*! \ingroup Channels
\brief The Additive White Gaussian Noise channel.

This class implements all functionalities of an Additive White Gaussian Noise channel.
User can set the \f$E_s/N_0\f$ (in dB or linear scale) and the seed used to initalize the random generator.
	
For an example of its use see e.g. the test program "test_LDPC.cpp".

\author Guido Montorsi
	*/
class AWGN_Channel  
{
public:
	//! Constructor
	AWGN_Channel(const int seedin=1  //!< Seed used to initalize the random generator.
		);

	virtual ~AWGN_Channel();

	//! Set the main parameters of block 
	void SetParameters(const double esn0,			//!< \f$E_s/N_0\f$ in linear scale
						const int seed=1			//!< Seed used to initalize the random generator. 
						);
	//! Run the AWGN channel
	void Run(const int ntics,						//!< Number of generated samples 	
			 const double* input,					//!< Input (real) signal 
			 double *out							//!< Output (real) signal
				);	
	
	//! Change the random generator seed
	void SetSeed(const int seedin			//!< Seed used to initalize the random generator
	){seed=seedin;};	

	//! Change Es/N0
	/*! Can be changed during simulation. The Es/N0 is in linear scale, 
	to set the EsN0 with dB quantities, see  Set_EsN0dB() */
	void Set_EsN0  (const double esn0				//!< \f$E_s/N_0\f$ in linear scale
);
	//! Change Es/N0
	/*! Can be changed during simulation. The Es/N0 is in dB, 
	to set the EsN0 with linear quantities, see  Set_EsN0() */
	void Set_EsN0dB(const double esn0			//!< \f$E_s/N_0\f$ in dB
	);	

	
#ifdef CTOPCOM
	void Run(const int ntics,						//!< Number of generated samples 	
		const cmplx* input,					//!< Input (real) signal 
		cmplx *out							//!< Output (real) signal
		)
	{
		Run(2*ntics,(const double*) input,(double*) out);
	}
#endif


	/** \cond INTERNAL */

	//! Run the AWGN channel on binary data with quantized LLR output  
	/*! Converts the input bits into antipodal values (+-1), add Gaussian noise and then converts
	the observed output into  LLR 
	\f[ \lambda=\frac{y_k}{\sigma_k}\f].
	LLR are then quantized.
	The output buffer can coincide with the input buffer */
	void Run(const int ntics,					//!< Number of generated samples 	
			 const int* input,					//!< Input signal 
			 int *out,							//!< Output signal 
			 double fact=8.						//!< Factor for quantization of LLR
				);

	//! Missing documentation
	void SetWrapped(const bool wrap=true){this->wrap=wrap;}


	//! Set the AWGN noise bias
	/*! The function pass a subset of index with biased value.
	This method is used for importance sampling techniques.
	*/
	void SetBias(
		const int K,			//! Length of bias subset
		const int* tset			//! Set of indexes
		);
	
	//! Position of bias subset
	void RunBias(
		const int ntics,	//! Number of generated samples
		const double*input,		//! Input signal
		double* out			//! Output signal
		);
	/** \endcond */

	int seed;		//!< Internal running seed used to generate the Gaussian distribution. 
	double level;	//!< Standard deviation of Gaussian random variables 

private:
	double gamma;
	double* computebasis(int K);
	int K;
	const int *tset;
	double* basis;
	double* nbuf;
	bool wrap;

};
