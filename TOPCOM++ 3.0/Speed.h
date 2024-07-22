// Speed.h: interface for the Speed class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdio.h>
#include <time.h>

/*! \file  
\brief Declaration of the class Speed
*/

/*! \ingroup SM
\brief Evaluates the speed of the simulation.

The measure is performed by counting the number of calls 
to the start() method per seconds.

The user should place a call to the start() method in the simulation loop.
Additionally, if desired, the user can put one more calls to the measure() method to evaluate
the amount of time spent on different sections of the loop

\author Guido Montorsi
*/ 
class Speed  
{
public:
	Speed();
	virtual ~Speed();

	//! Define the loop init point
	void start();

	//! Put a probing point.
	void measure();

	//! Display current measure.
	/* The routine measures the number of calls to the start()
	method for unit time. To display meaningful measures, specify the
	number of desired operations performed at each simulation step. e.g. number of bit
	processed...
	*/
	double Print(const double fact,	//!< Multiplying factor (units/loop)
		char* units,				//!< Name of units, e.g. "kbit/sec"
		FILE* stream=stdout			//!< output stream
		);
	void Reset();

private:
	clock_t ref;	//!< Reference time
	clock_t last;	//!< Last section call
	double *times;	//!< time spent in different sections 
	int ncycles;	//!< Number of cycles
	int n;
};
