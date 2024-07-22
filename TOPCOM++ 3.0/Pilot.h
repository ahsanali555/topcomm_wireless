// Pilot.h: interface for the Pilot class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of class Pilot
 */
/*! \ingroup Synchronization
\brief Insertion of pilot symbols into the data stream (see also Channel_Estimator).

The correspondent receiving block, that performs channel estimation based on the inserted
pilots is the Channel_Estimator.

The pilot pattern is specified through the method SetPilot as a
sequence of complex symbols to be inserted in the sequence of
modulated symbols (i.e. after the modulator but before filtering and
sampling). With this method the user can specify the period \f$P\f$ of
the pilot sequence and the complex value of the symbols in the
sequence \f$p\f$. These symbols may also be different from the
constellation points of the modulation.

The positions of the pilots in the data stream are instead specified through the method SetPositions(). 
In this method the user specifies
- The size \f$N\f$ of the fundamental data frame including the pilots
- The number \f$n\f$ of pilots inserted in the data frame
- The positions \f$k[i]\f$ where the pilots will be inserted (or added) in the frame. Obviously
    \f$k[i]\f$ must be in \f$[0,N-1]\f$. The vector \f$k[i]\f$ is assumed to be ordered.
Note that we let the possibility that the period \f$P\f$ of the pilot sequence  be different from the 
number \f$n\f$ of pilots
in a frame.

Two method of insertion are allowed. In the first, which corresponds to the traditional way of inserting pilots, 
the pilot symbol are inserted in the specified position "instead" of the data symbols. 
As a consequence the number of
data symbol in a frame is \f$N-n\f$ and the Pilot class increases the symbol rate by a factor \f$\frac{N}{N-n}\f$.

With the second insertion method the pilot symbol are inserted in the specified position "overlapped with"  
the data symbols. The data and pilot symbols are then simply added together in the specified positions. 
As a consequence
the number of data symbol in a frame is always \f$N\f$ and the Pilot class does not increase the symbol rate.

The method SetOverlapped allows to switch from one method of insertion to the other.

The method SetPower scales the amplitude of the pilot sequence (assumed to be originally with unit power), 
to give
the power specified by the user. The amplitude is computed also considering the positions 
where no pilots are inserted
and thus the scaling depends also on the ration \f$n/N\f$.

The Reset method allows  to reset the internal timing of the block and also to specify a 
possible initial offset
with respect to the frame boundary \f$N\f$. The default offset is 0 so that pilot insertion starts 
with a new frame and
with the first pilot in the pattern.

For an example of its use see e.g. the test program "test_iterative_phase.cpp".

\author Guido Montorsi
*/
class Pilot  
{
public:
	Pilot();
	virtual ~Pilot();
	//! Specify the positions of the pilots in the data stream 
	void SetPositions(const int N,	//!< Length of the reference data frame including the pilots
		const int npilots,			//!< Number of pilots inserted in the data frame
		const int* positions = 0		//!< Positions where the pilots will be inserted (or added) in the frame
		);

	//! Specify Regular structure 
	void SetRegularFrame(const int N,	//!< Number of  sections
		const int H,			//!< Length of Header 
		const int D,			//!< Length of Data section
		const int P			    //!< Number of Pilot Section
		);

	//! Specify the pilot pattern as a sequence of complex symbols to be inserted in the sequence of modulated symbols (default to end of frame)
	void SetPilot(const int period,		//!< Period of the pilot sequence 
		const double* pilot				//!< Complex values of the symbols in the sequence.
		);

	//! Insert the complex pilot sequence into a data stream
	int Run(int ntics,			//!< # of input complex samples (negative for setting the # of output samples)
		const double* data,		//!< Input data symbols
		double* datawithpilot	//!< Output data with pilot
		);

		//! Remove a complex pilot sequence from a data stream
	int RunInv(int ntics,				//!< # of input complex samples (negative for setting the # of output samples)
		const double* datawpilot,		//!< Input data symbols
		double* data					//!< Output data with pilot
		);

#ifdef CTOPCOM
	//! Insert the complex pilot sequence into a data stream
	//! Overload of Run with complex signals
	int Run(int ntics,			//!< # of input complex samples (negative for setting the # of output samples)
		const cmplx* data,		//!< Complex Input data symbols
		cmplx* datawithpilot	//!< complex Output data with pilot
		)
	{
	return (Run(ntics,(const double*) data,(double*)datawithpilot));
	}

	//! Remove a complex pilot sequence from a data stream
	//! Overload of Run with complex signals
	int RunInv(int ntics,				//!< # of input complex samples (negative for setting the # of output samples)
		const cmplx* datawpilot,		//!< Input data symbols
		cmplx* data						//!< Output data with pilot
		)
	{
	return (Run(ntics,(const double*) datawpilot,(double*)data));
	}
	
#endif

	//! Switch from insertion to overlapping of pilots and viceversa.
	void SetOverlapped(bool ov=true){overlap=ov;}

	//! Scale the amplitude of the pilot sequence (assumed to be originally with unit power)
	void SetPower(const double pow		//!< Power of the pilot sequence
		);

	//! Reset the internal timing of the block
	void Reset(int off=0		//!< Initial offset with respect to the frame boundary \f$N\f$
		);
	friend class Channel_Estimator;
	friend class Frame_Sync;
	friend class PLL;
	double fact;

	/** Return the vector of pilot of given length */
	double *GetPilotSequence(const int N);

private:
	int time;       // Internal clock
	int npp;	    // Counter of pilots positions (max npos)
	int npi;        // Counter of pilots inserted  (max pilots)
	int N;          // Periodicity of framing
	int npilots;    // Number of pilots in a frame
	int* positions;	// Positions of pilots
	int period;		// Period of pilot data
	double* pilot;	// Pilot data
	bool overlap;
};
