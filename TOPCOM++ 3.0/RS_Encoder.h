// RS_Encoder.h: interface for the RS_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

/*! \file 
\brief Declaration of class RS_Encoder.
*/

/*! \ingroup Codecs
\brief Reed-Solomon encoder

This class implements the functionalities of a systematic Reed-Solomon 
encoder. 
The parameters
of the RS code, that are introduced by the user, are:
- The codeword length \f$n\f$ in number of symbols(before shortening) 
- The information vector length \f$k\f$ in number of symbols (before shortening)
- The number of shortening symbols.

 The Reed-Solomon code generator polynomial and the GF generator polynomial 
can be also be specified by the user.
The GF generator polynomial is represented by an integer in hexadecimal format
(i.e. 0x11d represent the binary polynomial with coefficients 100011101:
\f$1+x^4+x^5+x^6+x^8\f$).
The code generator polynomial is in the form:
\f$g(x) = \prod_{j={\rm first\_root}}^{j<={\rm first\_root}+n-k}(x-a^{\rm fact}*j)\f$
i.e., the repeated product of the \f$n-k\f$ terms \f$(x-a^{\rm fact}*j)\f$ for values of \f$j\f$
from first_root to \f${\rm first\_root}+n-k\f$ inclusive. The
factor \f$a^{\rm fact}\f$ is a primitive root of GF\f$(2^m)\f$.

For an example of its use see e.g. the test program "test_RS.cpp".


\author Gabriella Bosco
*/

class RS_Encoder  
{
public:
	RS_Encoder();
	virtual ~RS_Encoder();

	//! Set the main parameters of the code
	void SetParameters(
			const int n,				//!< Codeword length (must be equal to \f$2^m-1\f$, with \f$m\in\{3,16\}\f$)
			const int k,				//!< Information word length (before shortening)
			const int nshort=0,			//!< Number of shortening symbols (default: 0)
			const int gfpoly=0,			//!< Generator polynomial of GF(2^m) (default: 0)
			const int first_root=1,		//!< Exponent of the first term in the product used in the construction of the code generator polynomial (default: 1)
			const int fact=1			//!< Factor used in the construction of the code generator polynomial (default: 1)
			);


	//! Run the encoder
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input symbols. 
		int* Output				//!< Output symbols.
		);
	
	
  
	friend class RS_Decoder;
	
	int kk;				//!< Information word length (w/o shortening).
	int ns;				//!< Number of shortening symbols
	int mm;             //!< Bits per symbol 
	int nn;             //!< Symbols per block (= (1<<mm)-1) 

private:
	void init_rs(int symsize,int gfpoly);
    int fcr;			//!< First root of RS code generator polynomial, index form
    int prim,iprim;		//!< Primitive element to generate polynomial roots
    int nroots;			//! <RS code generator polynomial degree (number of roots)
	int *index_of;      //!< Look-up table for the conversione of elements of \f$GF(2^m)\f$ from polynomial to index form
	int *alpha_to;		//!< Look-up table for the conversione of elements of \f$GF(2^m)\f$ from index to polynomial form
	int *genpoly;       //!< Generator polynomial 
	int *data,			//!< Information bits
		*bb;			//!< Redundancy bits
		
};
