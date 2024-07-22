#pragma once
#include <math.h>
/*! \file 
\brief Collects declarations of all basic C functions that are used to generate
random numbers
 */
/*! \ingroup functions
@{ */

//! Uniform generation of real  numbers  in ]0,1];
double	unif(
			 int&vc_un //!< Reference to the integer containing the seed.
						/*! Updated after each call */
			 );
//! Uniform generation of integers in ]0,2147483648];
int	unif_int(int  &vc_un//!< Reference to the integer containing the seed.
				 );
//! Box Muller transformation for Gaussian
inline double	Gaussian(int& seed  //!< Reference to the integer containing the seed.
				 );
//! Generate a random permutation.
int*	Random_Permutation(const int size, //!< Size of permutation.
						   const int seed  //!< Seed for the random generation.

);

inline unsigned int unif_int_qd(unsigned int &vc_un);
inline double   unif_qd(int &vc_un);


/*************************************************************************/
/****** Quick and dirty version *************/
/*************************************************************************/
inline double unif(int &idum)
{
	idum = 1664525*idum + 1013904223;
	return (unsigned int)idum/4294967295.;
}



//! Generate a Gaussian variable with unit variance 
inline double Gaussian(int& seed)
{
	double f1,f2,f3,a;
	static double x2=0.;
	static int isready=0;
	if(isready)
	{
		isready=0;
		return x2;
	}
	else
	{
		
		f3 = 100.0;
		f1=f2=0.;
		while (f3 > 1.0)
		{
			f1 = 2.0*unif(seed)-1.0;
			f2 = 2.0*unif(seed)-1.0;
			f3 = f1*f1 + f2*f2;
		}
		a = sqrt(-2.0 * log(f3)/f3);  
		
		x2 = f2 * a;
		isready=1;
		return (f1 * a);
	}
}


/*!@}*/

