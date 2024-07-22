// EyeDiagram.h: interface for the EyeDiagram class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
#include <stdio.h>
#include <stdlib.h>

/*! \file 
\brief Declaration of the class EyeDiagram
*/
/*! \ingroup Measurements 
\brief  Print the eye diagram of a signal.

The user can specifies the name of output file and the period of the waveform.

\author Guido Montorsi
*/


class EyeDiagram  
{
public:
	//! Constructor with parameters
	EyeDiagram(int ns,			//!< Period of waveform
		char* = NULL,		//!< Name of the output file
		bool iscmplx=false	//!< Is the input signal complex?
		);
	virtual ~EyeDiagram();
	
	//! Starts the printing of the eye diagram
	void Start(){start=true;ntime=0;};
	
	//! Ends the printing of the eye diagram
	void Stop(){fflush(file);printf("Eyediagram stored.\n");start=false;istart=-1;};
	
	//! Set starting and ending waveform
	void Start_Stop(
		const int startw,		//!< Starting instant 
		const int nwave=-1		//!< Number of periods
		);

	//! Runs the block
	void Run(const int, //!< Number of samples of the input signal
		const double*			//!< Input signal
		);
#ifdef CTOPCOM
	void Run(const int tics, //!< Number of samples of the input signal
		const cmplx*	input		//!< Input signal
		)
	{
		iscmplx=true;
		Run(tics, (const double*) input);
	}
#endif

private:
	bool iscmplx;		//!< Is the input signal complex?
	int ntime;			//!< Time index
	bool start;			//!< Has the storage started?
	int period;			//!< Period of waveform
	FILE* file;			//!< Name of the output file
	int	nwprint;		//!< Number of printed waveforms
	int istart;			//!< Starting instant 
	int nwave;			//!< Number of stored waveforms
};

