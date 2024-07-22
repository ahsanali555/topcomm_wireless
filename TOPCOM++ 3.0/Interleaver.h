#pragma once
#include <stdio.h>
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of class Interleaver and all its interface functions.
*/
/*! \ingroup Codecs
\brief Block and convolutional Interleaver (all types)

This class implements a generic interleaver, i.e., a
block that performs a permutation of the time axis of the input
sequence, providing at its output the same values in a different
order. The interleaver is described through a permutation \f$\pi\f$  of
size \f$N\f$ and can work with two different behaviors, block by block
or continuously. The two behaviors are implemented with the Run()
and RunCont() methods.

In the block-by-block behavior a set of \f$N\f$ input data is processed for
each tic. The output block is written by permuting the input block.
No memory is required and the latency is zero as a
block is produced for each input block.

In the continuous behavior a single data is processed for each tic. An
internal counter \f$i\f$ is incremented for each tic.
The data is written sequentially in an internal memory of size \f$N\f$
and one output is produced by reading the same memory in position \f$\pi(i)\f$.
At the beginning of the operation some
memory positions are read without having been written and this produces
a transient where some data are mixed with garbage initially stored in the memory.
All "convolutional" interleavers can be described by using such device.
This behavior implies an internal memory and
then requires synchronization.

Both methods can be called with a flag that specify if interleaving or de interleaving of data
is desired and can accept integer or real quantities.

Note that the back-to-back delay of a pair interleaver de-interleaver running in
continuous mode is always equal to \f$N\f$.

The "data" to be interleaved may also consist of a set of contiguous quantities,
giving the possibility of building
for examples "byte" interleavers.
The method SetCardinality() is used to set the number of quantities associated to
a single data, default assumes one data.

Several interface methods are available for this class that constructs permutations and interleavers for typical
applications.
  \see
	- Spread_Interleaver(),
	- Random_Interleaver(),
	- Congruential_Interleaver(),
	- UMTS_Interleaver(),
	- Interleaver_3GPP2(),
	- Row_by_Column(),
	- Forney(),
	- LoadFromFile().
	- Interleaver_InnDVBT()
	- Interleaver_SymbDVBT()

For an example of its use see e.g. the test program "test_interleaver.cpp".

\author Guido Montorsi
*/
class Interleaver
{
public:
	Interleaver();
	Interleaver(const int N, const int* perm=0);
	~Interleaver();

	//! Set the main parameters of the interleaver
	void SetParameters(
		const int N,		//!< Length of permutation
		const int *perm=0	//! Permutation (default to identity permutation)
		);


	//! Run the block interleaver with integer elements.
	void Run(const int ntics, //!< Number of interleaved blocks.
		const int *input,	//!< Input vector to be interleaved.
		int *output,		//!< Interleaved output vector.
		bool direct = true	//!< Direct or inverse interleaver.
		);

	//! Run the block interleaver with real  elements.
	void Run(const int ntics, //!< Number of interleaved blocks.
		const double *input,  //!< Input vector to be interleaved.
		double *output,			//!< Interleaved output vector.
		bool direct = true		//.!< Direct or inverse interleaver
		);



	//! Run the convolutional  interleaver with integer  elements.
	/*! The output buffer can coincide with the input buffer */
	void RunCont(
		const int ntics, //!< Number of interleaved blocks.
		const int *input,	//!< Input vector to be interleaved.
		int *output,		//!< Interleaved output vector.
		bool direct = true	//!< Direct or inverse interleaver
		);
	//! Run the convolutional interleaver with real elements.
	/*! The output buffer can coincide with the input buffer */
	void RunCont(
		const int ntics, //!< Number of interleaved blocks.
		const double *input,  //!< Input vector to be interleaved.
		double *output,			//!< Interleaved output vector.
		bool direct = true		//.!< Direct or inverse interleaver
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	//! Run the block interleaver with real  elements.
	void Run(const int ntics, //!< Number of interleaved blocks.
		const cmplx *input,  //!< Input vector to be interleaved.
		cmplx *output,			//!< Interleaved output vector.
		bool direct = true		//.!< Direct or inverse interleaver
		)
	{
	Run(ntics, //!< Number of interleaved blocks.
		(const double *)input,  //!< Input vector to be interleaved.
		(double *)output,			//!< Interleaved output vector.
		direct				//.!< Direct or inverse interleaver
		);
	return;
	}

		void RunCont(
		const int ntics, //!< Number of interleaved blocks.
		const cmplx *input,  //!< Input vector to be interleaved.
		cmplx *output,			//!< Interleaved output vector.
		bool direct = true		//.!< Direct or inverse interleaver
		)
		{
		RunCont(ntics, //!< Number of interleaved blocks.
		(const double *)input,  //!< Input vector to be interleaved.
		(double *)output,			//!< Interleaved output vector.
		 direct 		//.!< Direct or inverse interleaver
		);
		return;
		}

#endif


	//! Set the width of data (number of symbols).
	void SetCardinality(
		const int a	//!< Number of symbols in each data
		){card=a;}

	Interleaver* Inverse();

	void Display(FILE* =stdout);

	//! Save the permutation in a file
	void SaveToFile(char*);

	int  N;						//!< Length of permutation
	int  card;					//!< Cardinality of alphabet
	int* perm;					//!< Pointer to permutation.
	int time;					//!< Current value of internal clock

private:
	int* buffer;				//!< Buffer to store elements for convolutional interleaver.
	double* bufd;				//!< Buffer to store integer elements for convolutional interleaver.


friend
Interleaver* Spread_Interleaver(const int N,	
								int& S1,		
								int& S2,		
								const int& a,	
								int seed,		
								int trials,		
								bool linear		
								);


friend
Interleaver* Spread_Interleaver(const int N, 
								int seed	
								);	


friend
Interleaver* Random_Interleaver(const int N,	
								int seed	
								);

friend
Interleaver* Congruential_Interleaver(const int N, 	
									  int	
									  );

friend
Interleaver* UMTS_Interleaver(const int N		
							  );


friend
Interleaver* Interleaver_3GPP2(const int N		
							  );


friend
Interleaver* Interleaver_LTE(const int N		
							  );

friend
Interleaver* Interleaver_DVBSH(const int N		
							  );


friend
Interleaver* Interleaver_CCSDS(const int N		
							  );


friend
Interleaver* Interleaver_DVBRCS2(const int N		
								);


friend
Interleaver* Interleaver_InnDVBT(const int N		
							  );


friend
Interleaver* Interleaver_SymbDVBT(const int sh 
								  );


friend
Interleaver* Row_by_Column(const int M,					
						   const int N,				
						   const bool invcol		
						   );


friend
Interleaver* Row_by_Column(const int M,					
						   const int N,					
						   const int *colord			
						   );

friend
Interleaver* Forney(const int I,		
					const int N			
					);


friend
Interleaver* LoadFromFile(char*	name		
						  );

};


/*! \ingroup Interface
\brief Generate a spread random interleaver
*/

Interleaver* Spread_Interleaver(const int N,	//!< Interleaver size
								int& S1,		//!< Input spread
								int& S2,		//!< Output spread
								const int& a,	//!< Palindromia parameter
								int seed,		//!< Seed for random generation of permutations
								int trials,		//!< Number of attempts to perform before reducing the constraints
								bool linear		//!< Flag indicating if linear or square constraint are desired
								);

/*! \ingroup Interface
\brief Generate a spread random interleaver
*/

Interleaver* Spread_Interleaver(const int N, //! Interleaver size
								int seed=29861651	//!< Seed for random generation
								);	

/*! \ingroup Interface
\brief Generate an interleaver based on a ``random'' permutation of size \f$N\f$
*/

Interleaver* Random_Interleaver(const int N,	//!< Interleaver size
								int seed=76597	//!< Seed for random generation
								);
/*! \ingroup Interface
\brief Return an interleaver based on a permutation law of linear congruential type
*/

Interleaver* Congruential_Interleaver(const int N, 	//!< Interleaver size
									  int=57	//!< Seed for random generation
									  );
/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of UMTS standard
*/

Interleaver* UMTS_Interleaver(const int N		//!< Interleaver size
							  );

/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of the 3GPP2 standard
*/
Interleaver* Interleaver_3GPP2(const int N		//!< Interleaver size
							  );

/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of the 3GPP2 standard
*/

Interleaver* Interleaver_LTE(const int N		//!< Interleaver size
							  );
/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of the DVB-SH standard
*/

Interleaver* Interleaver_DVBSH(const int N		//!< Interleaver size
							  );

/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of the Deep space CCSDS  standard
*/

Interleaver* Interleaver_CCSDS(const int N		//!< Interleaver size (1784,3568,7136,8920)
							  );

/*! \ingroup Interface Standards
\brief Return an interleaver based on the permutation law of the DVBRC2 standard
*/

Interleaver* Interleaver_DVBRCS2(const int N		//!< Interleaver size (336,471,504,603,804,912,1200,1284,1539,1759,1884,2052,2259,3012)
								);

/*! \ingroup Interface Standards
\brief Return a bit interleaver as specified in DVB-T standard 300-744(inner interleaving)
*/

Interleaver* Interleaver_InnDVBT(const int N		//!< Interleaver size
							  );

/*! \ingroup Interface Standards
\brief Return the symbol interleaver specified in DVB-T standard 300-744 (inner interleaving)
*/

Interleaver* Interleaver_SymbDVBT(const int sh //!< sh=0 for long interleaver, sh=1 for short interleaver
								  );

 /*! \ingroup Interface 
\brief Return a row-by-column interleaver

Data written row by row and read column by column
*/

Interleaver* Row_by_Column(const int M,					//!< Number of rows of interleaver
						   const int N,					//!< Number of columns of interleaver
						   const bool invcol=false		//!< Invert rows order
						   );

 /*! \ingroup Interface
\brief Return a row-by-column interleaver

Data written row by row and read column by column
*/

Interleaver* Row_by_Column(const int M,					//!< Number of rows of interleaver
						   const int N,					//!< Number of columns of interleaver
						   const int *colord			//!< Ordering of columns 
						   );
/*! \ingroup Interface
\brief Returns a  Forney convolutional interleaver,
*/
Interleaver* Forney(const int I,		//!< Interleaving Number of branches \f$I\f$ 
					const int N			//!< Interleaving depth \f$N\f$
					);

/*! \ingroup Interface
\brief Load the permutation from a text file.

The user specifies the name of file. 
The first line of the files must be in the format "Size=%d\n" and specify the permutation size N.
Permutation is represented in following N lines with a number per line in the range [0,N-1]. The specified number represent the reading addresses
of a memory that is filled in the natural order.

*/
Interleaver* LoadFromFile(char*	name		//!< Name of file where permutation is loaded
						  );
Interleaver* PEG(const int N, int* partr, int *partl); // Sockets