// Fontan_Channel.h: interface for the Fontan_Channel class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include "Loo_Channel.h"
/*! \file 
\brief Declaration of the class Fontan_Channel
*/

/*! \ingroup Channels
\brief The Fontan Channel. Narrow band (single tap) Channel Model of Land Mobile Satellite.


  The channel is based on a three state model, with a Markov chain defining 
  the transition probabilities between the states. 
  Each state is characterized by a Loo_Channel with different parameters. 

  Parameters of the block are:
  - The 3x3 transition matrix describing the Markov chain
  - For each channel state,  the 3 parameters of the Loo model, namely the
  mean power of the lognormal LOS ([dB]), the standard deviation of the power of the lognormal
  LOS ([dB]), and the mean power of the diffusive rayleigh component
  - For each state, the minimal duration, measured in number of samples, of the state
  - The number of samples for the transitions between one state and the next. During transitions, 
  a linear combination of the two adjacent channels is performed
  - The doppler frequency, normalized to the sampling rate

  	For an example of its use see e.g. the test program "test_CompCodeComb_new.cpp"
	or its interface functions.

\author Guido Montorsi
  */
class Fontan_Channel  
{
public:
	Fontan_Channel();
	virtual ~Fontan_Channel();

	/*! Set the main parameters of the Fontanmodel */
	void SetParameters(
		const double *Pin,			//!< 3x3 Transition matrix of three state model
		const double *Pars,			//!< Parameters of Loo Model (alpha, Psi,MP) for the three channels
		const int* Nframe,			//!< Mimimal state length  (samples, for each state)
		const int Ncorr,			//!< Correlation length of shadowing process (samples)
		const int Ntrans,			//!< Transition length between states (samples)
		const double fd				//!< Normalized Doppler frequency (fd/Rs)
	);

	//! Generates samples of the Fontan process. 
	void Run(const int tics,  double* out);



	//! Generates samples of the Fontan process and multiplies them for an input signal.
	void	RunProd(const int tics, //!< Number of generated samples.
		const double* input,		//!< Input signal.
		double *output				//!< Output process.
		);

#ifdef CTOPCOM
	void Run(const int tics,  cmplx* out)
	{
		Run(tics,  (double*) out);
	}
	void	RunProd(const int tics, const cmplx* input,		cmplx *output)
	{
		RunProd(tics, 
			(const double*) input,	
			(double *)output				
			);
	}
#endif
private:
	Loo_Channel* Channel[3];
	double P[12];
	int Nframe[3];
	int currstate;
	int nextstate;
	int seed;
	int time;
	int Ntrans;





friend Fontan_Channel* Open_DVBSSP(const double Rs ,const double f0,const double speed); 
friend Fontan_Channel* Suburban_DVBSSP(const double Rs ,const double f0,const double speed); 
friend Fontan_Channel* Intermediate_Tree_Shadow_DVBSSP(const double Rs ,const double f0,const double speed); 
friend Fontan_Channel* Heavy_Tree_Shadow_DVBSSP(const double Rs ,const double f0,const double speed); 

};



/*! \ingroup Interface
\brief Generate a Fontan Channel as specified in TM-SSP0039r2, Open environment.
*/

Fontan_Channel* Open_DVBSSP(
				const double Rs ,	//!< Sampling rate of channel
				const double f0,	//!< Carrier frequency [GHz]
				const double speed	//!< Speed of mobile [km/h]
				);

/*! \ingroup Interface
\brief Generate a Fontan Channel as specified in TM-SSP0039r2, Suburban environment.
*/

Fontan_Channel* Suburban_DVBSSP(
				const double Rs ,	//!< Sampling rate of channel
				const double f0,	//!< Carrier frequency [GHz]
				const double speed	//!< Speed of mobile [km/h]
				);

/*! \ingroup Interface
\brief Generate a Fontan Channel as specified in TM-SSP0039r2, Intermediate tree shadow environment.
*/

Fontan_Channel* Intermediate_Tree_Shadow_DVBSSP(
				const double Rs ,	//!< Sampling rate of channel
				const double f0,	//!< Carrier frequency [GHz]
				const double speed	//!< Speed of mobile [km/h]
				);

/*! \ingroup Interface
\brief Generate a Fontan Channel as specified in TM-SSP0039r2, heavy tree shadow environment.
*/

Fontan_Channel* Heavy_Tree_Shadow_DVBSSP(
				const double Rs ,	//!< Sampling rate of channel
				const double f0,	//!< Carrier frequency [GHz]
				const double speed	//!< Speed of mobile [km/h]
				);


