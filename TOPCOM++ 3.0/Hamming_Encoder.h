// Hamming_Encoder.h: interface for the Hamming_Encoder class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
/*! \file 
\brief Declaration of class Hamming_Encoder.
*/
/*! \ingroup Codecs
\brief Extended Hamming encoder. 

  Extended Hamming codes are codes with correction capability one and a parity
  check for additional error detection. 

  For its use see the class PC_Encoder which has it as a member
  
\author Gabriella Bosco
*/

class Hamming_Encoder  
{
public:
	Hamming_Encoder();
	virtual ~Hamming_Encoder();

	//! Set the main parameters of the code
	void SetParameters(
		const int n1,      //!< Codeword length (must be equal to \f$2^m\f$)
		const int sh=0      //!< Shortened bits
		); 
	
	//! Run the encoder
	void Run(
		const int	blocks,		//!< Number of blocks.
		const int*	Input,		//!< Input signal. 
		int* Output			//!< Output signal.
		);

	friend class Hamming_Decoder;
	
	int n0,	//!< Codeword length.
		k0;	//!< Information word length.

	friend class PC_Decoder;

private:
	//!< Initialize the values of the coefficients of the generator and parity-check polynomials for Hamming code with lenght N=2^l-1.	
	void init(int l);
	//!< Divide a polynomial with coefficient in y by the generating polynomial g.	
	void genpoly_div(const int *y,const int n,int *q,int *r);
	//!< Create the generator polynomial for Hamming code with lenght N=2^l-1
	void polygen(int l);
	int sh;
	int *g;		//!< Generator polynomial
	int *h;		//!< Parity-check polynomial
	int *r;
	int *bufin,*bufout;
	int *tmp;			//!< Utility vector used in the routine genpoly_div
};

