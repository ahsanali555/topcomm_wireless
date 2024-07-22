// Long_Random_Generator.h: interface for the Long_Random_Generator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of class Long_Random_Generator.
*/
/*! \ingroup Random_Sources 
\brief Very long uniform random generator.

This random generator is based the algorithm "MWC1019" of Marsaglia, whose period exceeds \f$10^{9824}\f$.
MWC1019 has the property that every possible sequence of 1018 successive 32-bit integers will appear 
somewhere in the full period,  for those concerned with the "equi-distribution" 
in dimensions 2,3,...1016,1017,1018. The generator
requires 1019 initial 32 bit random values to start.

The user specifies the integer seed which is used to initialize the 1019 generators.
The two overloaded Run() methods provide both integers and doubles in the range [0,1].

For an example of its use see e.g. the test program "test_random_generators.cpp".

\author Guido Montorsi
 */
class Long_Random_Generator  
{
public:
	Long_Random_Generator();
	virtual ~Long_Random_Generator();
	
	//! Set the seed of the random generator.
	void SetParameters(const unsigned int seed  //!< Random generator seed
		);
	//! Generates a vector of integer random variables 
	void Run(const int ntics, //!< Number of generated random variables
				int* out	//!< Output vector
				);
	//! Generates a vector of double random variables in ]0,1] 
	void Run(const int ntics,   //!< Number of generated random variables 
			double* out		//!< Output vector
			);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	void Run(const int ntics,   //!< Number of generated random variables 
			cmplx* out			//!< complex Output vector
			)
	{
	Run(ntics,(double*) out);
	return; 
	}
#endif

private:
	unsigned long Q[1020];
	bool initialized;
	int p;
	static double RANDCONST;

};

