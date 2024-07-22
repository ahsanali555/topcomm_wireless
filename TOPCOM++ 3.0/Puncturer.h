// Puncturer.h: interface for the Puncturer class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of classes Puncturer and Depuncturer.
 */
/*! \ingroup Codecs
\brief Arbitrary puncturer or rate matching.

  The block allows to delete[] (puncture) a set of symbols from an input
  sequence according to some pattern that can be specified by the user.
  Alternatively, the method SetRateMatching allows to specify the number of
  input and surviving bits after puncturing. In this case a regular rate
  matching algorithm is applied that guarantees the desired puncturing ratio.

  The puncturing pattern is specified with a vector of integers. Each entry represent the number of
  times (possibly zero) that the considered input bit is replicated at the output. The puncturing 
  pattern is periodically applied to the input stream.

  For an example of its use see e.g. the test program "test_Convolutional.cpp".
\author Guido Montorsi
*/
class Puncturer  
{
public:
	Puncturer();
	virtual ~Puncturer();
	//! Set the puncturing pattern for arbitrary puncturing.
	void SetParameters(const int period,   //!< Length of puncturing pattern
		const int* pattern					//!< Puncturing pattern
		);

	//! Set the rate matching parameters (alternative to SetParameters()). 
	void SetRateMatching(const int ninp,  //!< Number of input symbols
		const int nout						//!< Number of output symbols
		);

	//! Run the puncturer on integer data (return number of generated outputs).
	int Run(const int tics,					//!< Number of processed inputs
			const int* input,				//!< Input symbols.
			int *out						//!< Output symbols.
			);

	//! Run the puncturer on double data (return number of generated outputs)
	int Run(const int tics,					 //!< Number of processed inputs
		const double* input,				//!< Input symbols
		double *out							//!< Output symbols
		);
#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,					 //!< Number of processed inputs
		const cmplx* input,				//!< Input symbols
		cmplx *out							//!< Output symbols
		)
	{
	return (Run(tics,(const double*) input,(double *)out));
	}
#endif


	//! Get the used period and puncturing pattern
	/*! Useful to get the puncturing pattern when rate matching algorithm is used.*/
	void GetPattern(
		int &period,		//!< Puncturing period (size of puncturing pattern)
		int* &pattern		//!< Puncturing Pattern
		);

	//! SISO (Bidirectional puncturer)
	/*! Returns the number of bits at the output of the puncturer. */
	int RunSISO(const int tics, //!< Number of processed input bits 
		const int* llri,   //!< llr on inp bits
		const int* llrc,   //!< llr on out bits		
		int* llruO,		   //!< ext on inp bits
		int* llrcO	   //!< ext on out bits.
		);
	int Nout(const int ninp);

	friend class Depuncturer;
	int period;
	int* pattern;
	int time;	//!< internal counter.
private:

	bool rm;	//!< Rate matching on.
	int num;	//!< Number of deleted/repeated bits.
	int den;	//!< Number of input bits.
};


//! Arbitrary de-puncturer (see also class Puncturer)
/*! \ingroup Codecs
  The block performs the inverse operation of Puncturer by inserting 0
  in the provided stream of integer in the position where the corresponding
  Puncturer has deleted a bit.

  Notice that the first input parameter specifies the number of generated outputs, so that
  its calling convention requires the same parameter of the corresponding Puncturer.

  For an example of its use see e.g. the test program "test_Convolutional.cpp".

\author Guido Montorsi
*/
class Depuncturer  
{
public:
	Depuncturer();
	virtual ~Depuncturer();

	//! Set the main parameters.
	void SetParameters(
		const Puncturer* pin			//!< Reference puncturer
		){p=pin;if(p->rm)time=1;else time=0;}

	//! Run the depuncturer on integer data (return number of processed inputs).
	int Run(
		const int tics,			//!< Number of generated outputs
		const int* input,		//!< Input symbols.
		int*out					//!< Output symbols.
		);

	//! Run the depuncturer on double data (return number of processed inputs).
	int Run(
		const int tics,			//!< Number of generated outputs.
		const double* input, 	//!< Input symbols.
		double*out				//!< Output symbols.
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,					 //!< Number of generated outputs
		const cmplx* input,				//!< Compelx Input symbols
		cmplx *out							//!< Complex Output symbols
		)
	{
	return (Run(tics,(const double*) input,(double *)out));
	}
#endif

	int time;
private:
	const Puncturer* p;
	


friend
Puncturer* Irregular_Puncturer(int K,				
							   int N,				
							   const int nr,		
							   const double* rates,	
							   const double* prob,	
							   int* len			
							   );


Puncturer* RC_Puncturer(int K,				
						int N				
						);
};

/*! 
Generate an irregular puncturing/repetition pattern with K input and N output.
The set of unnormalized rates and the correspondent unnormalized probabilities 
are passed through two vector. The vector are updated with the normalized values.
An additional optional vector is used to store the length of each puncturing
section.*/

Puncturer* Irregular_Puncturer(int K,				//!< Input bits
							   int N,				//!< Output bits.
							   const int nr,		//!< Number of rates
							   const double* rates,	//!< Rate of puncturing
							   const double* prob,	//!< relative length of puncturing section
							   int* len=0			//!< Length of puncturing section.
							   );

//! Rate Compatible puncturer 
/*! */
Puncturer* RC_Puncturer(int K,				//!< Input bits
						int N				//!< Output bits.
						);