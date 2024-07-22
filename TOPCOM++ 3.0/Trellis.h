// trellis.h: interface for the trellis class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <stdio.h>
/*! \file 
\brief Declaration of class Trellis.
*/
/*! \ingroup Codecs
\brief Description of a trellis.

The trellis is the basis structure for the description of convolutional codes.
The trellis is described through a set of
starting states \f$s\f$ and a set of input symbols \f$u\f$. 
A transition is characterized by the pair \f$(s,u)\f$. 
For each transition we associate the output symbol \f$c\f$ and the final state \f$s_f\f$ 
that must be in the same set of the starting state. 
The  available parameters are the number of states, the number of output symbols, the number of input symbols,
the number of input bits and output bits (only for binary encoders)


Two interface methods (Systematic_MaxRate() and Canonical()) are available for this class that 
construct trellises for typical applications.

  \see
  - Systematic_MaxRate()
  - Canonical()
  - Canonical2()
  - GetFromFile()

For an example of its use see e.g. the test program "test_Convolutional.cpp" 
or some other classes that use it as a member

\author Guido Montorsi
*/ 
class Trellis  
{
public:
	Trellis(const int=0,const int=0);
	Trellis& operator=(const Trellis& rhs);
//	void operator =(const Trellis & a);
	virtual ~Trellis();

	//! Return the time inverted trellis
	Trellis* Inverse();
	
	//! Return the terminating transition
	int Term(const int state	//!< Trellis state
		); 
	
    //! Return a vector that for each state gives the input symbol that drives the FSM toward the zero state. 
	/*! The method returns also the minimum number of trellis steps to go from any state to the zero state. */
	int* Termination(int& nterm  //!< Minimum number of termination steps
		) const;  

	//! Print the trellis structure to the desired output
	int DisplayTrellis(FILE* = stdout //!< Output stream
		) const;

	//! Computes the ending state cardinality
	int EndStateCard();

	//! Sort the trellis edges with respect to input labels
	void SortInput();



	//! Return a trellis with Hamming weight ofthe labels of this trellis as labels
	Trellis* Weight() const;

	//! Change the trellis encoder to make it systematic (if possible)
	bool SetSystematic();


	//! Returns the starting state associated to the pair (a= starting state ,b= input symbol)
	inline int starting(const int a,const int b)const {return trel[4*(a*ntrans+b)];} 

	//! Returns the input symbol associated to the pair (a= starting state , b= input symbol)
	inline int input(const int a,const int b) const {return trel[4*(a*ntrans+b)+1];} 

	//! Returns the output symbol associated to the pair (a= starting state , b= input symbol)
	inline int output(const int a,const int b)const {return trel[4*(a*ntrans+b)+2];}

	//! Returns the ending state associated to the pair (a= starting state ,b= input symbol)
	inline int ending(const int a,const int b)const{return trel[4*(a*ntrans+b)+3];}
	int* trel;		//!< Trellis description
	int nstate;		//!< Number of starting states
	int ntrans;		//!< Number of transitions for each state
	int no;			//!< Number of labels
	int k;			//!< Number of input bits
	int n;			//!< Number of output bits
//	int cons_len;	//!< Constraint length
	int* term;		//!< Termination vector
	int nterm;		//!< Number of terminating steps

friend Trellis* Systematic_MaxRate(const int k,const int* H);
friend Trellis* Canonical(int systematic,const int k,int n,const int* H,const int* Z);
friend Trellis* Canonical2(int systematic,const int k,int n,const int* H,const int* Z);
friend int ConstraintLength(const Trellis*);
friend Trellis* GetFromFile(char* name);
friend Trellis* Standard(const int a,const int b);



};


// Some common linear structures of convolutional encoders


//Trellis* LinearTrellis(const int,const int,const int,const int*,const int*,const int*,const int*);

/*! \ingroup Interface
\brief Return a systematic recursive encoder with \f$n=k+1\f$
*/

Trellis* Systematic_MaxRate(const int k, //!< Number of information bits
							const int* H //!< Feedback matrix
							);
/*! \ingroup Interface
\brief Return a canonical binary trellis encoder
*/

Trellis* Canonical(int systematic, //!< Is the trellis systemtic?
				   const int k,		//!< Number of information bits
				   int n,			//!< Number of encoded bits 
				   const int* H,	//!< Feedback matrix
				   const int* Z		//!< Feedforward matrix
				   );


/* \ingroup Interface
\brief Return a canonical binary trellis encoder
*/

Trellis* Canonical2(int systematic, //!< Is the trellis systemtic?
				   const int k,		//!< Number of information bits
				   int n,			//!< Number of encoded bits 
				   const int* H,	//!< Feedback matrix
				   const int* Z		//!< Feedforward matrix
				   );
/*\} */

//! Return the constraint length of a trellis 
int ConstraintLength(const Trellis*);


 /*! \ingroup Interface
\brief Return a trellis loaded from a file.
*/
Trellis* GetFromFile(char* name);

 /*! \ingroup Interface
\brief Return a trellis with standard delay shift structure.

The new state is obtained by left shift of the old state plus input label.
Output symbol is the concatenation of starting state and input symbol.
*/
Trellis* Standard(const int a,const int b);