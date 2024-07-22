// LDPC_Encoder.h: interface for the LDPC_Encoder class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <stdio.h>
/*! \file 
\brief Declaration of classes LDPC_Encoder and all its interface functions 
*/
/*! \ingroup Codecs*/
/*@{*/
/*!
\brief Low Density Parity Check Encoder 

The structure of the parity check matrix (in its sparse form)
is loaded from an external file, which has the following format:
- \f$n\f$ \f$m\f$
- \f$d_{c,max}\f$
- series of integer numbers defining the positions of the non-zero elements in row 1
- series of integer numbers defining the positions of the non-zero elements in row 2
- \f$\dots\f$
- series of integer numbers defining the positions of the non-zero elements in row \f$m\f$

where \f$n\f$ is the codeword length, \f$m=n-k\f$ is the number of rows of the parity check matrix,
\f$k\f$ is the length of the information vector, \f$d_{c,max}\f$ is the maximum number of non-zero entry in each row, 
corresponding to the maximum degree of check nodes.

The positions of the elements in the rows of the matrix are indicated by integer numbers in the range \f$(1,\;n)\f$.
If a row has less than \f$d_{c,max}\f$ elements, the corresponding line in the parity check matrix file is 
terminated by a zero.

The parity check matrix is assumed to be in lower triangular form, so that encoding operation can be performed
with linear complexity.

The interface function  LDPC_From_File(namefile) is used to generate an encoder for an LDPC code 
with the parity check matrix specified in the external file ""namefile".

The interface  function LDPC_DVBS2 generates an LDPC encoder 
implementing the DVB-S2 LDPC codes with block sizes 64Kbits and 16Kbits.

For an example of its use see e.g. the test program "test_LDPC.cpp".

\see LDPC_Encoder_CCSDS

\author Guido Montorsi
*/
class LDPC_Encoder  
{
public:
	LDPC_Encoder();
	virtual ~LDPC_Encoder();


	//! Construct the LDPC encoder
	void SetParameters(
		const int k,			//!< Number of input bits
		const int n,			//!< Number of output bits
		const int ncheck,		//!< Maximum degree of check nodes
		const int *H			//!< Parity check matrix (in lower triangular form)
		);

	//! Construct the LDPC encoder loading the structure from a file
	void SetParameters(const char* name);

	//! Reduce the PCM in a form suitable for linear time encoding
	/*! First makes the PCM in almost lower triangular form with gap bits.
	Then computes the dense matrix to generate the gap bits
	*/
	void MakeEncodable();

	//! Set the encoding algorithm (default to back substitution)
	/*!
		- 0: Back substitution, assume H in lower triangular form
		- 1: Use Generating matix
		- 2: WiMAX method  in 802.16e-2005  Matrix in H in circular block form
		- 3: Encoding with Compressed Generator matrix
		- 4: Method NG-LTE Qualcomm
		- 5: Gap  Matrix + Back substitution,
		*/
	void SetEncodingMethod(
		const int enctypein		//!< Encoding algorithm type
		){enctype=enctypein;}

	
	void Store_Matrix(const char * namefile) const;

	//! Optionally set a puncturing pattern
	void SetPuncturing(
		const int period,		//!< Period of puncturing pattern
		const int *pattern		//!< Puncturing pattern
		);

	// Optionally set a puncturing pattern
	void SetRegularRandom(
					const int K,			//!< Number of input symbols
					const int N,			//!< Number of output symbols (must be greater than K)
					double avdeg=0.			//!< Degree of variable nodes
					);


	void SetRegularRandom2(
					const int K,			//!< Number of input symbols
					const int N,			//!< Number of output symbols (must be greater than K)
					double avdeg=0.			//!< Degree of variable nodes
					);


	// Construct a random left irregular LDPC code
	void SetIrregularRandom(
		const int K,			//!< Number of input symbols
		const int N,			//!< Number of output symbols (must be greater than K)
		const int nterms,		//!< Number of terms in the left degree profile
		const double* ldeg,		//!< lefy degree values and relative percentage 
		bool edgeper=true		//!< Edge or node perspective	
		);

	// Construct a random left irregular LDPC code
	void SetLeftIrregularRandom(
		const int K,			//!< Number of input symbols or 
		const int N,			//!< Number of output symbols (must be greater than K)
		const int ndv,			//!< Number of terms in the left degree profile
		int* dv,			//!< left degree values 
		const double* pdv,		//!< left degree percentage 
		const int Z = 1,			//!< Quasi-Cyclic extension?	
		const bool edgeper = true,  	//!< Edge or node perspective	
		const bool PEG = false		//!< PEG optimization?	
		);

	// Construct a random right irregular LDPC code
	void SetRightIrregularRandom(
			const int K,
			const int N,
			const double avdeg,
			const int ndc,
			const int* dc,
			const double* pdc,
			bool edgeper=true
				 );

	//! Run the encoder using the selected encoding algorithm (default to back substitution)
	/*! Several encoding methods are implemented depending on the structure of the parity check matrix.
	The default one assumes that parity check matrix is in lower triangular form, so that back substitution can be
	performed.
	*/
	virtual void Run(const int ntics,	//!< Number of blocks
				const int* inp, //!< Input bits
				int* out		//! Output bits
				);

	//! Print the parity check matrix
	void Display(
		FILE* file=stdout,   //!< Output stream
		const int Z=0        //!< Assume Quasi cyclic structure with lifting Z
		) const;

	void GetCode(int & N, int & E, int & Z, int & dC, int *& P, int *& add, int *& sh, int *& deg) const;


	//! Check if the Parity check matrix is lower triangular
	/*! Returns a boolean indicating if the parity check matrix is lower triangula and consequently can
	be encoded with back substitution */
	bool IsLowerTriangular() const;

	//! Return a flag indicating if the provided codeword is valid
	bool IsCodeword(const int* inp) const;

	//! Shorten the code
	/*! A set of input bits are  set to zero, and the 
	same bits are removed from the generated codeword.
	Previous value is superseeded
	*/
	void Shorten(const int Xs, //!< Number of shortened bits					
		const int last=0
	)
	{
		this->lastsh = last;
		K   = K1 - Xs+shorten;
		N  -= (Xs-shorten);
		shorten=Xs;
	}

	// Activate the HARQ, data are stored for future retransmission
	void SetHARQ();

	void SetVarPerm(const int type=0);

	int * PEG(const int nedges, const int K, const int * degv, const int * degc);
	void Repeat(const int repeat,	//!< Repetition number
				int* out					//! Output bits
				);

	void LLRCombine(const int repeat,	
				int* llr,					
				const int* out				
				);
	void AddDegreeOneVN(const int Ndegone);
	int SetMultiEdge(const int N, const double * NU,const double * MU, const int Z=1);

	int K;		//!< Number if input bits to the encoder
	int N;		//!< Output bits (after possible puncturing)

	int K1;		//!< Number of input bits (including shortened bits)
	int N1;		//!< Length of parity check matrix
	int period;	//!< Period of puncturing pattern
	int* pattern; //!< Punturing pattern

	//! friends interfaces
	friend LDPC_Encoder* LDPC_From_Generating_Matrix(char *,char *);
	friend LDPC_Encoder* LDPC_From_File(char *);

	friend LDPC_Encoder* LDPC_DVBS2(const int , const int );
	friend LDPC_Encoder* LDPC_DVBSX(const int, const int);
	friend LDPC_Encoder* LDPC_WIMAX(const int rate, const int N, const bool A);
	friend LDPC_Encoder* LDPC_WiFi(const int rate, const int N);
	friend LDPC_Encoder* LDPC_WiFi80211ad(const int rate, const int length, const bool lift, int* seed);
	friend LDPC_Encoder* LDPC_RegularRandom2(const int K,const int N,	int ncheck);
	friend LDPC_Encoder* LDPC_NRQUALCOMM(const int K,const int N, int& family,const int sol);
	friend LDPC_Encoder * LDPC_NR5G(const int K, const int N, int &family, int &cluster);
	friend LDPC_Encoder * LDPC_NR5G_strong(const int K, const int N, int &family, int &cluster);
	friend int * Stronger(const int k, const int n, const int * H, const int st1, const int st2);
	friend class CLDPC_Encoder;
	friend class CLDPC_Decoder;
	friend class LDPC_Decoder;
	friend class MODM_LDPC_Decoder;

	void RunM(const int ntics,const int* inp, int* out);
	int M;		//!< Cardinality of alphabet (default to 2)

	int GetZ() const { return z; }

	// Computes the matrix to generate p2 bits
	int z;			//!< Expansion factor (when available)
	int enctype;	//!< Type of encoding algorithm
	int shorten;	//! Flag to indicate if shortening is applied

private:
	void ComputeGapMatrix(int g,const int z=1);
	int * PEG(const int K, const int N, const int nedges, const int *degv, const int* degc);

	// Some buffers
	int* buff;		//!< Internal buffer for puncturing
	int* Hconn;		//!< Connections of parity check bits
	int* Gconn;		//!< Compressed generator matix
	int* G;			//!< Generating matrix
	int* p;			//!< Permutation of input
	int* bufHARQ;	//!< buffer for storing coded data
	int* gg;		//!< Matrix for generating gap bits (Richardson-Algorithm)
	int* order;		//!< degree order for selective HARQ

	int gap;        //!< Gap bits (Richardson-Algorithm)
	int lastsh;     //! Flag for shortening at the end of inf block (default to zero);
	int ncheck;		//!< Maximum degree of check nodes
	int nmaxvar;	//!< Maximum degree of variable nodes (for compressed generator matrix)
	int x;			//!< circular shift of v(0) (for LTE encoder)
	int hhh;		//!< For Qualcomm LTE encoding
	int aaa;		//!< For Qualcomm LTE encoding
	bool HARQ;
	int counter;
	int ccore;      //! For NR LTE codes, core of the PCM
	int* permvar;   //! optional additional permutation on variable nodes
};

/**************************** Interfaces **********************/
/*! \ingroup Interface 
\brief Build an LDPC starting from the parity check matrix stored in  a file

See the manual for an explanation of conventions
\author Guido Montorsi
*/
LDPC_Encoder* LDPC_From_File(
							 char * namefile //!< Name of the file
							 );

/*! \ingroup Interface 
\brief Build an LDPC starting from the generating and parity check matrix.

  The systematic part must be at the end of the parity check matrix
\author Guido Montorsi
*/
LDPC_Encoder* LDPC_From_Generating_Matrix(
							 char * namegenmat, //!< Name of the file storing the generating matrix
							 char * nameparmat    //!< Name of the file storing thhe parity check matrix
							 );

/*! \ingroup Interface Standards 
\brief Return an LDPC as specified in the WiMAX 802.16e-2005 standard.
  The encoding algorithm is that specified in the appendix of the standard (method 1a).
  \author Guido Montorsi
*/
LDPC_Encoder* LDPC_WIMAX(
						 const int rate,	//!< Rate of encoder (12,23,34,56)
						 const int N,		//!< Code length in bits. Must be a multiple of 96 in [576,-]
						 const bool A=true  //!< A code ? (only for rate 23 and 34)
						 );

/*! \ingroup Interface Standards
\brief Return an LDPC as specified in the WiFi 802.11-2012 standard.

  \author Guido Montorsi
*/
LDPC_Encoder* LDPC_WiFi(
						 const int rate,	//!< Rate of encoder (12,23,34,56)
						 const int N		//!< Code length in bits. Must be a multiple of 96 in [576,-]
						 );

/*! \ingroup Interface 
\brief Return an LDPC as specified in the WiFi 802.11ad-2012 standard.

  \author Guido Montorsi Standards
*/
LDPC_Encoder* LDPC_WiFi80211ad(
	const int rate,				//!< Rate of encoder (12,58,34,1316,78(not standard))
	const int length=672,		//!< Length of codeword (multiple of 672 are allowed but not in the standard)
	const bool lift = false,		//!< Samsung lifting for generating length 1,344
	int *seed=0	//!< seed for rando lifting generation
	);

/*! \ingroup Interface Standards
\brief Return an LDPC as specified in the DVB-S2 standard.
\author Guido Montorsi
*/
LDPC_Encoder* LDPC_DVBS2(
						/*! Rate of the code
						 - 14 : 1/4 (0,1,2)
						 - 13 : 1/3 (0,1)
						 - 25 : 2/5 (0,1)
						 - 12 : 1/2 (0,1,2)
						 - 35 : 3/5 (0,1)
						 - 23 : 2/3 (0,1)
						 - 34 : 3/4 (0,1,2)
						 - 45 : 4/5 (0,1)
						 - 56 : 5/6 (0,1)
						 - 89 : 8/9 (0,1)
						 - 910: 9/10(0)
						 */
						 const int type,	
						 const int shortblock=0 //!< short code (0:64,800 1:16,200 2:4096)
						 );

/*! \ingroup Interface Standards
\brief Return an LDPC as specified in the new set of DVB-SX standard.
\author Guido Montorsi
*/
LDPC_Encoder* LDPC_DVBSX(
								/*! Rate of code
								long:
								15,29,1345,920,1120,
								2645,2845,2336,2536,1318,
								79,90180,96180,100180,104180,
								116180,124180,128180,132180,135180,
								140180,144180,150180,154180,1830,2030,2230

								short: 1145,415,1445,715,815,2645,3245
								*/

								const int rate,
								const int shortblock//!< short code (0:64,800 1:16,200)
								);

int* Gaussian_Elimination(unsigned int *H, const int n, const int k);

LDPC_Encoder* LDPC_RegularRandom2(const int K,const int N,	int ncheck=0);

LDPC_Encoder* LDPC_NRQUALCOMM(
	const int K,		//!< Information block size[bits]
	const int N,		//!< Codeword size [bits]
	int &family,		//!< Required code family (1,2,3) default to optimum family according to the rate
	const int sol=-1	//!< Required code solution, default to that minimizing Shortened and Punctured bits
	);
/*@}*/
