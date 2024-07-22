// Transponder.h: interface for the Transponder class.
//
//////////////////////////////////////////////////////////////////////

#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#pragma once
#include <stdio.h>
#include "Filter.h"
#include "NonLinearity.h"
/*! \file
\brief Declaration of the class Transponder and all its interfaces
*/

/*! \ingroup Channels
\brief Satellite Transponder, including IMUX, NonLinearity and OMUX.

This class embeds three objects of type Filter, NonLinearity and Filter and simulates the behavior 
of a complete Satellite transponder.
The objects are initialized by the user and their pointer are passed through the method SetParameters 
to configure the block.
A pointer to null corresponds to the absence of the correspondent block.
Notice that the detructor deletes the embedded blocks although they are not allocated by the constructor.
The Reset() method reset the contents of IMUX or OMUX filters if they are present. 

\see
 - EDMT_Transponder()
 
 \author Guido Montorsi
*/
class Transponder  
{
public:
	Transponder();

	// Destructor
	/* The destructor deletes the embedded blocks (IMUX, OMUX and NonLinearity)*/
	virtual ~Transponder();

	//! Configure the Transponder by passing the pointers to embedded sub blocks
	void SetParameters(
		Filter* IMUX,			//!< Pointer to IMUX filter (set to zero if absent)
		NonLinearity* TWT,		//!<  Pointer to NonLinearity (set to zero if absent)
		Filter* OMUX			//!<  Pointer to OMUX filter (set to zero if absent)
		);	
	
	//! Run the Transponder for a given number of complex samples
	/*! The output buffer can coincide with the input buffer */
	void Run(const int ntics,			//!< Number of complex input samples
		const double *inp,				//!< Input signal (complex)
		double*out						//!< Output signal (complex)
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
		void Run(const int ntics,			//!< Number of complex input samples
		const cmplx *inp,				//!< Input signal (complex)
		cmplx* out						//!< Output signal (complex)
		)
		{
		Run(ntics,(const double*)inp,(double*)out);
		return;
		}
#endif


	//! Reset the memory contents of IMUX and OMUX filters (if present)
	void Reset()
	{
		if(withIMUX)IMUX->Reset();
		if(withOMUX)OMUX->Reset();
	}


	Filter* IMUX;			//!< Pointer to embedded IMUX filter
	Filter* OMUX;			//!< Pointer to embedded OMUX filter
	NonLinearity* TWT;		//!< Pointer to embedded nonlinearity
	bool withIMUX;			//!< Flag to activate the IMUX 
	bool withOMUX;			//!< Flag to activate the OMUX 
	bool withTWT;			//!< Flag to activate the Transponder

	//! Pointer to total impulse response of transponder (zero if not available)
	/*! Convolution of IMUX and OMUX  impulse responses */
	double* hresp;			

	//! number of samples of impulse reponse
	int nresp;
	double Powout,ap;
	double Powinp;


private:


/*! \ingroup Interface
\brief Generate a transponder as specified in the SOW of ESA contract EDMT

This Interface is used to create an instance of class Transponder with features described in
	- "Enhanced Digital Modem Techniques Development and Validation. STATEMENT OF WORK" Ref. TEC-ETC/2010.75/NA
	Issue: 1, Revision 0 30.05.2010

The three embedded blocks can be activated independently with the correspondent flags.

The 3dB bandwidth of the transponder can be tuned with the parameter Bw. The total impulse 
response  of the transponder, neglecting non linearity, is computed first.Depending on the active blocks, this impulse response can correspond to
the IMUX only, the OMUX only, or both filters.

The time axis is then properly scaled by resampling the impulse responses of filters to achieve the
desired bandwidth. The resulting total impulse response is stored in the public buffer #hresp.  


Four types of nonlinearities are possible:
		- 0 Ideal_NonLinearity()
		- 1 DVBS2_Linearized()
		- 2 DVBS2_NotLinearized()
		- 3 TWT_EDMT()
*/
friend
Transponder* EDMT_Transponder(
		bool withIMUX,		//!< Flag to activate IMUX
		bool withOMUX,		//!< Flag to activate OMUX
		bool withTWT,		//!< Flag to activate OMUX
		int NLtype,			//!< Type of nonlinearity (0,1,2,3)
		double Bw			//!< 3dB bandwidth of transponder (normalized to sampling rate)
		);

friend
Transponder* DVBSX_Transponder(
		bool withIMUX,		//!< Flag to activate IMUX
		bool withOMUX,		//!< Flag to activate OMUX
		bool withTWT,		//!< Flag to activate OMUX
		int NLtype,			//!< Type of nonlinearity (0,1,2,3)
		double Bw			//!< 3dB bandwidth of transponder (normalized to sampling rate)
		);

};

Transponder* EDMT_Transponder(bool withIMUX,bool withOMUX,bool withTWT,int NLtype,double Bw);
Transponder* DVBSX_Transponder(bool withIMUX,bool withOMUX,bool withTWT,int NLtype,double Bw);