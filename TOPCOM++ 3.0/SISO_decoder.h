// SISO_Decoder.h: interface for the SISO_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Trellis.h"
#include "Delay.h"
#include "Include.h"
/*! \file 
\brief Declaration of class SISO_Decoder.

*/
/*! \ingroup Codecs 
\brief Soft Input Soft Output decoder for time invariant trellis.

This general block is the core of several iterative decoding blocks that 
use convolutional codes. It computes output soft values (extrinsic information) in the form of LLR 
of input and output sequences of a trellis encoder starting from the corresponding sequence of
soft values. For more details on the algorithm see the manual.

For an example of its use see e.g. the test program "test_CPM.cpp" or some other decoder classes
that use it as a member (PCCC_Decoder, SCCC_Decoder...)


\author Guido Montorsi
*/	
class SISO_Decoder  
{
public:
	SISO_Decoder();
	virtual ~SISO_Decoder();

	//! Set the main parameters of the SISO block
	void SetParameters(const Trellis*,
		const int nini	 = 100 , //!< Training window size 
		const int ngroup = 100,	 //!< Updating window size 
		const double f	 = 8.	 //!< Factor for representation of LLR.
		);

	//! Run the SISO over a block of data, assuming termination.
	int Run(
		const int nsteps,		//!< Number of trellis steps.
		const int* Inf,			//!< LLRs on information bits. (Set to zero if not available).
		const int* Cod,			//!< LLRs on coded bits.	(Set to zero if not available).
		int* OInf,				//!< EXTs on information bits (Set to zero if not required).
		int* OCod				//!< EXTs on coded bits		  (Set to zero if not required).
		);

	//! Run the SISO continuously over a block of data (sliding window).
	void RunCont(
		const int nsteps,	//!< Number of trellis steps.
		const int* Inf,		//!< LLRs on information bits. (Set to zero if not available).
		const int* Cod,		//!< EXTs on information bits (Set to zero if not required).
		int* OInf,			//!< EXTs on information bits (Set to zero if not required).
		int* OCod,			//!< EXTs on coded bits		  (Set to zero if not required).
		const int* index=0	//!< Optional index of tranmitted symbol for M.I. evaluation
		);


	//! Set the precision factor for quantized input.
	void SetPrecision(
		const double f //!< Precision factor.
		);

	//! Set the input LLRs to be on binary symbols
	void SetBinaryInput();

	//! Set the output LLRs to be binary on binary symbols
	void SetBinaryOutput();

	//! Set the input and output LLRs to be on minimal alphabets.
	/** The input and output alphabet cardinality is decomposed into primes and input
	and output alphabets are assumed to be defined on sub-alphabets with 
	the correspondent cardinalities.
	*/
	void SetMaxSplit();

	//! Set the  decomposition of the input alphabet to the sub-alphabets specified by the user.
	void SetInpSplit(
		const int,			//!< Number of sub-alphabets.
		const int* InpSplit	//!< Cardinality of sub-alphabets.
		);
	//! Set the  decomposition of the output alphabet to the sub-alphabets specified by the user.
	void SetOutSplit(const int,	//!< Number of sub-alphabets.
		const int* CodSplit		//!< Cardinalities of sub-alphabets.
	);

	//! Get the current mutual information rate
	double Get_MI();
	//! Reset MI meter
	void Reset_MI();

	int *ic; //!< Look-up table for max* operator
	bool APPc;//!< Flag for returning TOTAL information instead of EXTRINSIC information
	bool APPu;//!< Flag for returning TOTAL information instead of EXTRINSIC information

	//!< Reset the decoder
	void Reset();
void OneStep(
	   const int *A, const int* Inf,const int* Cod,const int *B,  
	   int *OA, int* OInf,  int* OCod,int *OB);

	//! Final state assumed for backward recursion
	int finalstate;
	int Py,Pyx;
	int tsteps;
	bool extnorm;
	bool resetMI;

	bool finalterm;


private:

	/*! \brief Run a single window of BCJR with given initialization of backward recursion. 
		 Only information on coded symbols
		 */
	void Window(
		const int*, 
		int*, 
		int*,
		int=0); 

		/*! \brief Run a single window of BCJR with given initialization of backward recursion.
	Keep previous values for forward recursion
		 Information on both coded symbols and information symbols
	*/
	void Window(
		const int*,		
		const int*, 
		int*, 
		int*,
		int=0);

	inline	int maxx(const int* buffer,const int);
	int splittedI;       // Flag for splitted input
	int splittedO;       // Flag for splitted output
	int nsi;
	int nso;
	int nlfi;
	int nlfo;
	int* SplitI;
	int* SplitO;
	int* tempbufferI;
	int* tempbufferO;
	double factor;
	
	void Collect(
		const int* Input,	// Input LLF
		int* Output,		// Output LLF
		const int ns,		// Number of symbols
		const int* split,	// Cardinality of each alphabet
		const int nllf,		// Total number of input LLF
		const int off=0		// Offset for input circular buffers
		);

	int* Project(int nsym, int ns,const int* split, const int nlf, const int* LLRI, int* LFO, int* LLRO,int*);
	void Ext_Comp(const int nsym,  const int ns,  const int *split,const int nlf, int*& inp, int*& out);

	int* backward;
	int* inputback;
	int* outputback;

	int* forward;
	int* inputforw;
	int* outputforw;

	int* infpartition;
	int* codpartition;

	int* beta_buffer;
	int* beta_swap;
	int* alfa_swap;
	int* maxx_buffer;
	int* buff_o;
	int* buff_i;

	int* apath;
	int* napath;

	int* gammaU;
	int* gammaC;
	int* gammaUo;
	int* gammaCo;

	int maxtab;

	int ns;     
	int ni;
	int no;
	int nini;
	int ngroup;
	int nti;
	int nto;
	int time;
	int *buffi,*buffo;
	Delay* del;
};