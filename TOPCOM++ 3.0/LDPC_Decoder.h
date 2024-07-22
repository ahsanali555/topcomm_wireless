// LDPC_Decoder1.h: interface for the LDPC_Decoder class.
//
//////////////////////////////////////////////////////////////////////


#pragma once

#include "LDPC_Encoder.h"
/*! \file 
\brief Declaration of class LDPC_Decoder
*/
/*! \ingroup Codecs
\brief Low Density Parity Check Decoder based on belief propagation.

The user specifies:
- The maximum number of iterations \f$I_{max}\f$ in the decoding algorithm.  
- A factor \f$f\f$ (the same used in the Demodulator class), which defines the
 correspondance between 
actual LLRs values \f$\lambda_i\f$ (double) and their quantized version \f$\lambda_Q\f$ (integer):
\f[
\lambda_Q=\lfloor\lambda_i\cdot f+0.5\rfloor.
\f]

  The Run() method provides optimal decoding based on belief propagation.
  The faster RunSubopt() method performs decoding using a suboptimal implementation
  of the check node operation (minimum with offset on magnitude). 

For an example of its use see e.g. the test program "test_LDPC.cpp".

\author Gabriella Bosco, Guido Montorsi
*/
  
class LDPC_Decoder  
{
public:
	LDPC_Decoder();
	virtual ~LDPC_Decoder();

	//! Set the main parameters of the block 
	/*! The code parameters are specified through the corresponding encoder */
	void SetParameters(const LDPC_Encoder*, //!< Pointer to the encoder
		const int niter,				    //!< Number of iterations
		const double fact=8.				//!< Precision factor
		);


	//! Set the decoder to work with reinforced BP 
	/*!  */
	void SetRBP(
		const double* gamma		//!< Reinforcing parameter vector (one for each iteration)
		);


	//! Set the decoder to work with "Damping" 
	/*!  */
	void SetDamping(
		const double damp		//!< Reinforcing parameter vector (one for each iteration)
		);


	//! Run the decoder (optimal and slower version)
	double Run(const int ntics,				//!< Number fo blocks 
		const int* LLRh,					//!< Input LLRs of coded bits
		  int* out,							//!< Decoded LLRs of data bits.
		  const int *data = 0,				//!< Data bits for genie-aided stopping rule.
		  int* out2=0,						//!< Optional additional buffer to store output codeword LLRs
		  bool reset=true					//!< Flag to reset the extrinsic information at beginning of operations
			);

	//! Run the sub-optimal (min-sum-offset) decoder
	double RunSubopt(const int ntics,				//!< Number fo blocks 
		const int* LLRh,					//!< Input LLRs of coded bits
		  int* out,							//!< Decoded LLRs of data bits.
		  const int *data = 0,				//!< Data bits for genie-aided stopping rule.
		  const int beta =3,				//!< offset value (def to Optimal value for a factor of 8 in LLR)
		  int* out2=0,						//!< Optional additional buffer to store output codeword LLRs
		  bool reset=true					//!< Flag to reset the extrinsic information at beginning of operations
			);

	//! Run the sub-optimal (damped) decoder
	double RunDamped(const int ntics,		//!< Number fo decoded codewords 
		const int* LLRh,					//!< Input LLRs of coded bits
		  int* out,							//!< Decoded LLRs of data bits.
		  const int *data = 0,				//!< Data bits for genie-aided stopping rule.
		  const int beta =3,				//!< offset value (def to Optimal value for a factor of 8 in LLR)
		  int* out2=0,						//!< Optional additional buffer to store output codeword LLRs
		  bool reset=true					//!< Flag to reset the extrinsic information at beginning of operations
			);

	//! Run the LDPC decoder with layered scheduling
	double RunLayered(const int ntics,		//!< Number fo decoded codewords 
		const int* LLRh,					//!< Input LLRs of coded bits
		  int* out,							//!< Decoded LLRs of information bits.
		  const int *data = 0,				//!< Optional pointer to Data bits for genie-aided stopping rule.
		  bool opt=true,					//!< Flag for optimal or sub-optimal (min-sum-off) decoder
		  const int beta =3,				//!< Offset value (def to Optimal value for a factor of 8 in LLR)
		  int* out2    =0,					//!< Optional additional buffer to store output codeword LLRs
		  bool reset=true					//!< Flag to reset the extrinsic information at beginning of operations
			);


	//! Displays the variable and check degree distribution (node perspective) of the LDPC code.
	void Display(FILE* file = stdout			//!< Output file (default to stdout)
	) const;



	//! Displays the variable and check node degrees of the LDPC code.
	void LongDisplay(
		FILE * file = stdout //!< Output file (default to stdout)
	) const;
	
	//! Run decoder with non binary input LLR
	double RunLayeredQ(const int ntics, const int * LLRCh, int * out, const int * data, bool opt, const int beta);

	int maxx(const int a, const int b);


	//! Return 
	//! Set the precision of llrs representation
	void Set_Precision(const double fact);

	double GetDensity() const { return (double)this->nones / N1; }
	int sat;		//!< Saturation of messages (Optional) 
	int niter;			//!< Number of iterations
	int nones;			//!< Number of ones in the parity check matrix
	bool stop;			//!< Flag to activate stopping rule;
	int checkmax;       //!< Verified checks for stopping rule 
	friend double GetThreshold(const LDPC_Decoder* dec, const double S, const int P);

private:
	int* nconncheck;	//!< Variable node distribution
	int* nconnvar;		//!< Check node distribution
	const LDPC_Encoder* encoder;
	int K;				//!< Input bits
//	int N;				//!< Output bits
	int N1;				//!< Variable bits
	int* perm;			//!< Permutation that define the connections
	int maxdeg;			//!< Maximum degree of check nodes
	int* buf;			//!< buffer for check nodes operation
	int* EXT;		    //!< internal buffer
	inline int g(const int a, const int b);	//!< g operator
	int*  ic;	//!< LUT for  max*
	int maxtab; //!< size of LUT
	double damp;
	const double* gamma;		//< Reinforcing parameter
	int* OLD;
	int* Lambda;

};
