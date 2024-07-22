// NonLinearity.h: interface for the NonLinearity class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdio.h>
#include "Modulator.h"
/*! \file
\brief Declaration of the class NonLinearity and all its interfaces
*/

/*! \ingroup Channels
\brief Generic memoryless non linearity, specified through AM/AM and AM/PM curves.


SetParameters()  lets the user specify the required
AM/AM and AM/PM curves through a set of points. The user specifies:
- The number of provided points in the AM/AM and AM/PM curves
- The value (in dB) of the input power at specified points
- The value (in dB) of the output power at specified points
- The value (in degrees) of the output phase rotation at specified points

The set of points are then interpolated with splines in order to build a new 
table with very high resolution. The
look-up table contains, in step of 0.01dB of the input power, 
the AM/AM and AM/PM characteristics, from the minimum to
the maximum specified powers.

In the Run method the power (in dB) of the instantaneous input power is evaluated and quantized to the closest
value available in the table. The output amplitude and phase is then efficiently computed using by a LUT access. 
If the input power fall outside the LUT range, linear extrapolation is performed using the 
points on the boundaries.

The method SetIBO() reduces the input power by the specified
amount of dBs. The provided value would then coincide with the true
IBO only if the original signal power is one.

Interfaces
\see
 - DVBS2_Linearized()
 - DVBS2_NotLinearized()
 - XBand_SSPA()
 - XBand_TWT()
 - SSPA_MHOMS()
 - TWT_MHOMS()
 - TWT_ESA()

 See "test_nonlinearity.cpp" or its interface functions for an example of its use. 

 \author Guido Montorsi
*/
class NonLinearity  
{
public:
	NonLinearity();
	virtual ~NonLinearity();

	//! Specify the required AM/AM and AM/PM curves
	void SetParameters(
		const int npAM,			//!< Number of provided points in AM/AM curve
		const double *inpam,	//!< Value (in dB) of the input power at specified points (AM curve)
		const double *outam,	//!< Value (in dB) of the output power at specified points
		const int npPM,			//!< Number of provided points in  AM/PM curves
		const double *inppm,	//!< Value (in dB) of the input power at specified points for (PM curve)
		const double *outpm,		//!< Value (in degrees) of the output phase rotation at specified points
		const double stepdb=0.01	//!< Step in dB of the internal LUT for AM/AM and AM/PM curve
		);	
	
	//! Run the nonlinarity
	/*! The output buffer can coincide with the input buffer */
	void Run(const int ntics,			//!< Number of complex input samples
		const double *inp,				//!< Input signal (complex)
		double*out						//!< Output signal (complex)
		);

	//! Reduces the input amplitude by the specified amount of dBs
	void SetIBO(
		const double IBOin,		//!< Input Back off in dB
		bool phacom=true		//!< Flag to compensate phase at working IBO
		);

	//! Scale input and output by the specified factor
	void SetGain(const double Gin){G=Gin;}

	//! Print on the desired output the AM/AM and AM/PM curves as computed by the interpolation algorithm.
	void Display(FILE* file=stdout);

	//friend NonLinearity* Inverse(const NonLinearity* ref);
	//friend NonLinearity* Ideal_NonLinearity();

#ifdef CTOPCOM
	void Run(const int ntics,			//!< Number of complex input samples
		const cmplx *inp,				//!< Input signal (complex)
		cmplx*out						//!< Output signal (complex)
		)
	{
		Run( ntics,(const double *)inp,(double*)out);
	}
#endif
	friend class EGRET_Channel;
private:
	double G;
	double IBO;
	double *amsin;
	double *amcos;
	double *amam;
	double min;
	double max;
	double delta;
	int nsamp;

	/* Extrapolation coefficients */
	double	gAM, oAM, gPM,oPM;



/*! \ingroup Interface
\brief Generate the linearized amplifier specified in the standard DVB-S2
*/
friend
NonLinearity* DVBS2_Linearized();

/*! \ingroup Interface
\brief Generate the nonlinarity corresponding to the not linearized 
amplifier specified in the standard DVB-S2
*/
friend
NonLinearity* DVBS2_NotLinearized();


/*! \ingroup Interface
\brief Generate the Xband SSPA  specified by ESA.
*/
friend
NonLinearity* XBand_SSPA();

/*! \ingroup Interface
\brief Generate the Xband TWT  specified by ESA.
*/
friend
NonLinearity* XBand_TWT();

/*! \ingroup Interface
\brief Generate a solid state power amplifier. 
  The AM/AM and AM/PM curves were taken from SSPA model of the MHOMS project.
*/
friend
NonLinearity* SSPA_MHOMS();

/*! \ingroup Interface
\brief Generate a traveling wave tube amplifier
  The AM/AM and AM/PM curves were taken from the TWTA model MHOMS project.
*/
friend
NonLinearity* TWT_MHOMS();

/*! \ingroup Interface
\brief Generate a TWT specified by ESA.
*/
friend
NonLinearity* TWT_ESA();


/*! \ingroup Interface
\brief Generate a TWT specified by ESA for a carrier frequency of 26GHz.
*/
friend
NonLinearity* TWT_ESA26GHz();



friend
NonLinearity* Ideal_NonLinearity();

friend
NonLinearity* Inverse(const NonLinearity* ref);

/*! \ingroup Interface
\brief Generate a TWT specified in the ESA-EDMT contract.
*/
friend
NonLinearity* TWT_EDMT();

/*! \ingroup Interface
\brief Generate the ground HPA power amplifier specified in DVB-SX Channel Models
*/
friend 
NonLinearity* GroundHPA_DVBS2();

/*! \ingroup Interface
\brief Generate a Solid State Power amplifier for KU band.
*/
friend
NonLinearity* KuBand_SSPA();

friend
NonLinearity* TWT_Rapp(const double p,const double a);

/*! \ingroup Interface
\brief Generate a static predistorter, which impose a predefined output constellation.
*/
friend
Modulator* StaticPredistorter(const Modulator*, NonLinearity *TWT);


};

NonLinearity* DVBS2_Linearized();

NonLinearity* DVBS2_NotLinearized();

NonLinearity* XBand_SSPA();

NonLinearity* XBand_TWT();

NonLinearity* SSPA_MHOMS();

NonLinearity* TWT_MHOMS();

NonLinearity* TWT_ESA();

NonLinearity* TWT_ESA26GHz();

NonLinearity* Ideal_NonLinearity();


NonLinearity* Inverse(const NonLinearity* ref);


NonLinearity* TWT_EDMT();
NonLinearity* GroundHPA_DVBS2();
NonLinearity* KuBand_SSPA();
NonLinearity* TWT_Rapp(const double p,const double a);

Modulator* StaticPredistorter(const Modulator*, NonLinearity *TWT);
