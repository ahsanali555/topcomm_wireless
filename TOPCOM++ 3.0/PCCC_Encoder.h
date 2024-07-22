// PCCC_Encoder.h: interface for the PCCC_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "Trellis.h"
#include "Interleaver.h"
/*! \file
\brief Declarations of PCCC_Encoder class and related interface functions
 */
/*! \ingroup Codecs 
\brief Generic binary Parallel Concatenated Convolutional Code encoder.

The user specifies the trellises of the two encoders and the interleaver. The
trellises must have input and output symbols belonging to the sets \f$U\f$ and \f$C\f$ with cardinalities
that are powers of two \f$U=2^k\f$ and \f$C=2^n\f$, so that we
can associate to each trellis step a number of input bits and output bits. 
Typical trellises are then obtained using binary encoders as those constructed 
through the interface method Canonical.

The PCCC's constituent encoders work in terminated mode and the interleaver 
is a block interleaver so that the
PCCC as a whole is a block encoder.


Optionally, the user can add a puncturing block to both encoders through the methods SetPunctUpper and
SetPunctLower. The puncturing are specified through their periodicity and the puncturing pattern 
as a vector of
"1" and "0", where "1" means bit not punctured and "0" means bit punctured.
In order to apply puncturing the bits generated from the encoder are taken from top 
to bottom on a single trellis
section and then sequentially. 
The pattern is applied to all encoded bits but the terminating ones.

For an example of its use see e.g. the test program "test_PCCC.cpp".

  \see Other Interfaces: 
  - PCCC_UMTS_Encoder()
  - PCCC_3GPP2_Encoder()
  - PCCC_Encoder_CCSDS()

\author Guido Montorsi
*/
class PCCC_Encoder  
{
public:
	PCCC_Encoder();
	virtual ~PCCC_Encoder();
	//! Set the parameters of teh decoder
	void SetParameters(const Trellis* upper,	//!< Trellis upper encoder.
		const Trellis*lower,					//!< Trellis lowrer encoder.
		const Interleaver*  interleaver //!< Interleaver
		);				
	//! Add a puncturing block to the upper encoder
	void SetPunctUpper(const int per, //!< Puncturing period
		const int* pat					//!< Puncturing pattern
		);
	//! Add a puncturing block to the lower encoder
	void SetPunctLower(const int per, //!< Puncturing period
		const int* pat				//!< Puncturing pattern
		);

	//! Special method added to deal with unusual terminations techniques like 3GPP2
	/*! See PCCC_3GPP2_Encoder() */
	void SetPunTerm(
		const int peru,				//!< Puncturing period termination upper.
		const int* patu,			//!< Puncturing pattern upper.
		const int perl,				//!< Puncturing period termination lower.
		const int* patl				//!< Puncturing pattern lower
		);

	
	//! Run the encoder (parallel output)
	/*! A tic of the Run method accepts a block of \f$K\f$ information bits and produces 
	a block of \f$N\f$ coded bits. All encoded bits of the upper encoder are stored in the first part of the 
	output buffer, while all encoded bits of the lower encoder are stored in the final part */
	void Run(const int tics,					//!< Number of encoded blocks
		const int* inp,							//!< Input bits
		int* out								//!< Output bits
		);

	void RunBGAN(const int tics, const int * input, int * out);

	void SetSequential();

	int K;			//!< Input bits
	int Nu;			//!< Output bits upper encoder
	int Nl;			//!< Output bits lower encoder
//	int N;			//!< Output bits

	friend class PCCC_Decoder;
	friend class TX_BGAN_102744;

private:
	void CodeSize();
	int punperu;
	int punperl;

	int *order; //!< Change order of output bits
	int pertermu;
	const int*  pattermu;
	int perterml;
	const int*  patterml;

	const int* patternu;
	const int* patternl;
	const Trellis *upper;
	const Trellis *lower;
	const int* perm;
	int ntu;		//!< Number of trellis steps upper encoder
	int ku;			//!< Number of input bits  upper encoder
	int nu;			//!< Number of output bits upper encoder
	int ntl;		//!< Number of trellis steps lower encoder
	int kl;			//!< Number of input bits  lower encoder
	int nl;			//!< Number of output bits lower encoder
	int *tu;			//!< input symbol for the termination of upper encoder
	int *tl;			//!< input symbol for the termination of lower encoder
	int ntermu;			//!< Number of steps for the termination of upper encoder
	int nterml;			//!< Number of steps for the termination of lower encoder

	int punctBGAN;	//!< Set a special puncturing according to the BGAN standard


friend
PCCC_Encoder* PCCC_UMTS_Encoder(
	const int K			
						);


friend
PCCC_Encoder* PCCC_3GPP2_Encoder(
	 const int rate,	
	 const int K		
		);





friend
PCCC_Encoder* PCCC_Encoder_CCSDS(
	 const int rate,	
	 const int K		
								 );


friend
PCCC_Encoder* PCCC_DVBSH_Encoder(
	 const int rate,	
	 const bool data,		
	 const bool comp);

};



/*! \ingroup Interface
Generate a PCCC encoder as specified in the UMTS standard.
*/

PCCC_Encoder* PCCC_UMTS_Encoder(
	const int K			//!< Information blocks size (39<K<5114)
						);

/*! \ingroup Interface 
Generate a PCCC encoder as specified in the 3GPP2 standard.
*/

PCCC_Encoder* PCCC_3GPP2_Encoder(
	 const int rate,	//!< Rate of encoder. (e.g. 12=1/2) available 34,23,35,12,13,14,15 
	 const int K		//!< Information block size. (127<K<32768).
		);




/*! \ingroup Interface 
Generate a PCCC encoder as specified in the CCSDS standard (102..
*/

PCCC_Encoder* PCCC_Encoder_CCSDS(
	 const int rate,	//!< Rate of encoder. (e.g. 12=1/2) available 12,13,14,16 
	 const int K		//!< Information block size. (1784,3568,7136,8920).
								 );

/*! \ingroup Interface 
Generate a PCCC encoder as specified in the DVB-SH standard (ETSI EN 302.83).
*/

PCCC_Encoder* PCCC_DVBSH_Encoder(
	 const int rate,	//!< Rate of encoder. (e.g. 12=1/2) available 15,29,14,27,13,25,12,23}
	 const bool data=true,		//!< Flag for data (K=12282) or signal (K=1146).
	 const bool comp=false		//!< Flag for complementary code.
		);

/** \cond INTERNAL */
/*! \ingroup Interface 
Generate a PCCC encoder as specified in the 3GPP (Rel.9 LTE) standard.
*/
//friend
//PCCC_Encoder* PCCC_LTE_Encoder(const int K);

/** \endcond */