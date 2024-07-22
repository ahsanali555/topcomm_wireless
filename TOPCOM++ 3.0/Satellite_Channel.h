#pragma once
#include "Delay.h"
#include "ShiftFrequency.h"
#include "Filter.h"
#include "Transponder.h"
#include "Phase_Jitter.h"
#include "FIFO.h"

#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of the class Satellite_Channel
*/

/*! \ingroup Systems Channels
	

\brief Satellite Channel Model as specified in the document "DVB‐SX Channel Models", by 
the TM-S2 Channel Model Group, with Chairman and editor: A. Ginesi. April 2014. 

The user is assumed to be familiar with this document.
The channel model includes the following set of blocks and impairments

- Non linear gateway amplifier
- Uplink Interference from adjacent carrier gateways
- Satellite: IMUX, Satellite on board amplifier, and OMUX
- Downlink cochannel interference
- Downlink adjacent carriers interference
- Fase and frequency offsets
- Amplitude distortions


All channel blocks modelling the different impairments can be inserted or removed toggling the correspondent flags  with the SetParameters() method.


\author Guido Montorsi
*/
class Satellite_Channel
{
public:
	Satellite_Channel(void);
	~Satellite_Channel(void);

	//! Configure the satellite channel.
	/*! Each flag corresponds to the activation of one impairment in the model.
	The number of input samples equals the number of output samples.*/
	void SetParameters(
		const bool ampup,			//!< [in] Gateway HPA
		const bool intup,			//!< [in] Interference on uplink
		const bool imux,			//!< [in] IMUX
		const bool twt,				//!< [in] Transponder amplifier
		const bool omux,			//!< [in] OMUX
		const bool intdown,			//!< [in] Downlink interference
		const bool phasenoise,		//!< [in] Phase Noise
		const bool freqoff,			//!< [in] Frequency offset
		const bool ampldist,		//!< [in] Amplitude distortions at RX
		const int ns=6,				//!< [in] Samples per symbol
		const double Rs  =34,		//!< [in] Baud rate [MHz]
		const double psi = 36./34.,	//!< [in] Ratio of Tranponder bandwidth to baud rate
		const double fT  = 40./34.,	//!< [in] Ratio of carrier spacing to baud rate
		const double APU =0.,		//!< [in] Gain of reference channel w.r.t. adjacent channels for uplink interference [dB]
		const bool staggered=true,	//!< [in] Staggered downlink cochannel interference
		const bool deepfading=false,//!< [in] Flag for deep fading
		const int phmask=1			//!< [in] Type of phase noise mask
		);

	// Run the satellite simulation
	void Run(
		const int ntics,	  //!< [in] Number of samples
		const double* input,  //!< [in]  Channel input samples(complex baseband)
		double* out			  //!< [out] Channel output samples (complex baseband)
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(
		const int ntics,	  //!< [in] Number of samples
		const cmplx* input,  //!<  [in] Channel input samples(complex baseband)
		cmplx* out			  //!< [out] Channel output samples (complex baseband)
		)
	{
		Run(ntics,(const double*) input,(double*) out);
		return;
	}
#endif
	ShiftFrequency* fs[9];
	Delay* delup;
	Delay* deldown;
	Delay* delamp;
	Transponder* TTT;
	NonLinearity* HPA;
	Phase_Jitter* ph;
	Phase_Jitter* phlow;
	AtoD* Upsamp;
	FIFO* fifo;
	Filter* cable;

	double* phstore; //!< Pointer to an external buffer to store the generated phase samples.	
	double Powout;	 //!< Running average Power at the output of TWT (after OMUX)
	double ap;		 //!< Parameter of the one pole filter to compute the running average power Satellite_Channel::Powout

private:
	int ns;			//!< Number of samples per symbol
	double psi;		//!< Transponder bandwidth to baud rate ratio
	double gamma;	//!< Frequency spacing to baud rate ratio
	double Ap;
	bool staggered;
	bool deepfading;
	double Df;
	double tau;
	double ripple;
	double delta;
	double Rs;		//!< Baud rate
	double apu;		// Gain of reference channel w.r.t. the adjacent channels (linear)
	void SetAPU(
		const double APU	//!< Gain of reference channel w.r.t. adjacent channels [dB]
		)
	{
		this->apu=pow(10.,APU/20.);
	}
	void SetPhaseNoise(
		const int PNtype //!< Phase noise mask type
		);

	bool intup;
	bool ampup;
	bool imux;
	bool twt;
	bool omux;
	bool intdown;
	bool phasenoise;
	bool freqoff;
	bool ampdist;
	double *buff1;
	double *buff2;
	double *buff3;
	double *phase,*phase1,*phase2;
	static const int nmax;


};

