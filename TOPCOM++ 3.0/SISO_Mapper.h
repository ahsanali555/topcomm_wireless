// SISO_Mapper.h: interface for the SISO_Mapper class.
//
//////////////////////////////////////////////////////////////////////


#pragma once

/*! \file 
\brief Declaration of class SISO_Mapper and its interface functions.
*/ 
/*! \ingroup Codecs
\brief Generic Soft Input Soft Output Module for arbitrary mappings.


This general block is the core of several iterative decoding blocks. 
It computes output soft values (extrinsic information) in the form of LLR 
of input and output sequences of a generic memoryless mapping, starting from 
the corresponding sequence of soft values. 

The mapping is a correspondence not necessarily invertible between an input set that is a
cartesian product of some sub-symbols and an output set which is also  a
cartesian product of some sub-symbols.

Typical mapping is that existing between the bits and the constellation points of a modulator.

For more details on the algorithm see the manual.

For an example of its use see e.g. the test program "test_iterative_phase.cpp" 
or some other classes that use it as a member (Demodulator...)

\author Guido Montorsi
 */

class SISO_Mapper  
{
public:
	SISO_Mapper();
	virtual ~SISO_Mapper();

	//! Set the parameters of the SISO mapper
	void	SetParameters(
		const int *mapping, //!< Input/Output Mapping 
		const int ni,		//!< Number of input alphabets
		const int *ipart,	//!< Cardinality of input alphabets
		const int no, 		//!< Number of output alphabets
		const int *opart,	//!< Cardinality of input alphabets
		const double fact=8.//!< Factor for LUT computation of max* 
		);

	//! Set the precision factor
	void SetPrec(const double fact);
	
	//! Run the SISO mapper with sequential inputs.
	/** This method allows to pass the input/output LLRs as a vector of pointers to the vectors
	containing the LLRs of constituent subalphabets.
	Each pointer is associated to one input or output alphabet. If one pointer is set to zero it means
	that the LLRs of the correspondent input/output sub-alphabet is not provided/required.
	*/
	void Run(
		const int tics,			//!< Number of symbols
		const int **luI,		//!< Vector of pointers to the input LLRs  on information symbols
		const int **lcI,		//!< Vector of pointers to the input LLRs  on coded symbols
		int **luO,				//!< Vector of pointers to the output LLRs  on information symbols
		int **lcO				//!< Vector of pointers to the output LLRs  on coded symbols
		);

	//! Run the SISO mapper with parallel inputs.
	/** This method allows to pass the input/output LLRs as pointers to the vectors
	containing the LLRs of ALL subalphabets. The LLRs of the sub-alphabets are assumed to 
	be stored sequentially in a single vector.
	If one pointer is set to zero it means that 
	the LLRs of ALL the input/output sub-alphabets are not provided/required.*/
	void	RunParall(
		const int tics,			//!< Number of symbols
		const int *luI,			//!< Pointer to the input LLRs  on all information symbols
		const int *lcI,			//!< Pointer to the input LLRs  on all coded symbols
		int *luO,				//!< Pointer to the output LLRs  on all information symbols
		int *lcO				//!< Pointer to the output LLRs  on all coded symbols
		);
	
	int  nli;		 //!< Total Number of inputs llr
	int  nlo;		 //!< Total Number of outputs llr	

	friend class Iterative_Decoder;

/* Interface functons */

/*! \ingroup interface
Returns a SISO mapper for the "natural" mapping.

The mapping is computed as follows. 
A prime factor decomposition is performed over the integer N:
\f[ N=\prod_{j=1}^m p_j\f]

A generic element \f$n \in N\f$ is represented as follows
\f[ n=\sum_{j=1}^m n_j \prod_{i=1}^{j-1}p_j  \f]
which  identifies the  following mapping
\f[ n \leftrightarrow (n_1,\ldots,n_m) \f].
between a set of m sets with cardinalities equal to the 
factors of N, and the set of integers N. 
In the particular case where N is a power of 2, the mapping coincides with the classical
binary representation of integers. 
*/
friend SISO_Mapper* Natural(
		const int N, //!< Cardinality
		const bool direct,		//!< Specifies if bit to symbol (true) o symbol to bit mapping is considered
		const double fact		//!< Factor for LUT computation of max* 
		);

/*! \ingroup interface
Returns a SISO mapper for a "Gray" mapping.
\sa Natural()
*/
	friend SISO_Mapper* Gray(
		const int N,				//!< Cardinality
		const bool direct,		//!<  Specifies if bit to symbol (true) or symbol to bit mapping is considered
		const double fact		//!< Factor for LUT computation of max* 
		);

private:
	int* mapping;   // Mapping
	int* inpcar;	 //!< Cardinality of each input alphabet
	int* outcar;	 //!< Cardinality of each output alphabet
	int  ni;		 //!< Number of inputs
	int  no;		 //!< Number of outputs	
	int  Ni;		 //!< Total cardinality

	int* inp;
	int* out;
	int* normu;
	int* normc;
	inline void maxx(int &a, const int b);
	int* ic;
	int maxtab;

};
SISO_Mapper* Natural(
	const int N, //!< Cardinality
	const bool direct=true,		//!< Specifies if bit to symbol (true) o symbol to bit mapping is considered
	const double fact=8.		//!< Factor for LUT computation of max* 
	);

SISO_Mapper* Gray(
	const int N,				//!< Cardinality
	const bool direct=true,		//!<  Specifies if bit to symbol (true) or symbol to bit mapping is considered
	const double fact=8.		//!< Factor for LUT computation of max* 
	);
