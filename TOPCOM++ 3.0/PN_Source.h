// PN_Source.h: interface for the PN_Source class.
//
//////////////////////////////////////////////////////////////////////

#pragma once


#include "Include.h"
/*! \file 
\brief Declaration of the class PN_Source
*/
/*! \ingroup Random_Sources 
\author Guido Montorsi
\brief Pseudo casual binary source based on Maximal length sequences. 
 */
class PN_Source  
{
public:

	//! Set the parameters of the PN source.
	PN_Source(unsigned int=31,  //!< log-2 of the length of generated PN sequence (degree of p.polyinomial).
		unsigned int seed=1		//!< Initial content of shift register.
		);
	virtual ~PN_Source();

	//! Set a polynomial defined by the user.  
	void SetFeedback(
		const uint64 fb			//!< Feedback polynomial coefficients (LSB is the lowest power coefficient)
		);


	//! Change the contents of the shift register
	void SetSeed(unsigned int seedin	//!< New seed.
		){seed= (uint64)(seedin&mask);}

	//! Run the PN source generating random bits
	void Run(
		const int tics,			//!< Number of generated bits
		int* output				//!< Pointer to buffer to store generated bits
		); 

	
	//! Run the inverse PN source  
	void RunInv(
		const int tics,			//!< Number of generated bits
		int* output				//!< Pointer to buffer to store generated bits
		); 

	//! Run the PN source and scramble data
	void Scramb(
		int tics,				//!< Number of scrambled bits
		const int* input,		//!< Pointer to input bits
		int* output,			//!< Pointer to scrambled bits
		bool inv =false,		//!< Flag for inverse source
		bool soft=false			//!< Flag for soft inputs
		); 

	__int64 seed;		//!< Shift register current content

private:
	unsigned int     degree;
	uint64 feedback;
	uint64 mask;
	static const uint64 poly[64];
};
