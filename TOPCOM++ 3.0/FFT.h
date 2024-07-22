// FFT.h: interface for the FFT class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file 
\brief Declaration of the class FFT
*/
/*! \ingroup DSP 
\brief Performs Fast Fourier Transform on complex sequences.

For an example of its use see e.g. the method "Spectrum::Run()"

\author Guido Montorsi
*/
class FFT  
{
public:
	//! Constructor
	FFT(const int n=12		//!< Base-2 logarithm of the maximum size of FFT performed.
		);
	virtual ~FFT();

	//! Run the on-place FFT on a complex signal stored in a single vector
	/*! Overloaded to accept complex type */
	void Run(const int m,	//!< Base-2 logarithm of the size of FFT performed.
		double *io,			//!< Buffer of complex data to be transformed on place.  
		const double s=1	//!< Select direct (1) or inverse (-1) transform 
		);

	//! Run the on-place FFT on a real signal stored in a single vector
	void RunR(const int m,	//!< Base-2 logarithm of the size of FFT performed.
		double *io,			//!< Buffer of complex data to be transformed on place.  
		const double s=1	//!< Select direct (1) or inverse (-1) transform 
		);

	//! Run the on-place FFT on a complex signal represented as two vectors containing real and imaginary parts
		void Run(const int m,	//!< Base-2 logarithm of the size of FFT performed.
		double *real,			//!< Buffer of real part of data to be transformed on place.  
		double *imag,			//!< Buffer of imaginary part of data to be transformed on place. 
		const double s=1		//!< Select direct (1) or inverse (-1) transform
		);
	static double* W;
	static int* inver;

#ifdef CTOPCOM
	void Run(const int m,	//!< Base-2 logarithm of the size of FFT performed.
		cmplx *io,			//!< Buffer of complex data to be transformed on place.  
		const double s=1	//!< Select direct (1) or inverse (-1) transform 
		)
	{
		Run(m,	
		(double *)io,		
		s	
		);
	}

#endif
private:
	static int nmax;	// Maximum number of bits
	static int NMAX;
	static void ComputeTables(const int N);
};



