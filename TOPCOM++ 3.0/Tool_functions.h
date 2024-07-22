/*! \file 
\brief Collects declarations of some basic functions and algorithms that are used in TOPCOM++
*/

/*! \ingroup functions 
@{ */

//TOOLS.H
#pragma once

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Include.h"

int gcd(int a, int b);
int mcm(int a, int b);

//! Greatest common divisor
void extended_gcd(int a, int b, int& lastx, int& lasty);

//! Computes Hamming weight of (the binary representation) an integer
inline int Hamming(unsigned int a)
{
	int i=0;
	while(a>0)
	{
		i+=(a&1);
		a>>=1;
	}
	return i;
}

	const double PI=3.1415926535897932384626433832795;

#ifdef WIN32


	//! Complementary error function
	//double	erfc(double xx);

	//! Set a the priority of the process to "low".
	void LowPriority();

	//! Check heap consistency
	int heapdump( void );
#endif

//! Inverse of the complementary error function
double	inverfc(const double xx);

//! Return the binomial of a on b
double	binom(const double a, const double b);
double	binom2(const double a, const double b);
double	binom3(const double a, const double b);


//! Prime factor decomposition of N
int*	prime(const int N);

//! Return the log of the binomial of a on b
double	logbinom(const double a,const double b);
//! Return the log of the binomial of a on b
double	logbinom2(const double a,const double b);
//! Return the log of the binomial of a on b
double	logbinom3(const double a,const double b);
//! Return the max* of a and b
///int		maxx(const int a, const int b);

//! Perform log interpolation between two numbers
double	log_interpolate(double ystar, double* y, double* x, int n);
//! Perform  interpolation between two numbers
/*! Returs the absissa yielding the required value of y */
double interpolate(double ystar,  const double* x, const double* y,const int n);

int		log2(int a);
int		getfromfile_i(FILE* file,char* name);
float	getfromfile_f(FILE* file,char* name);
double	getfromfile_d(FILE* file,char* name);

//! MAX* transform
void MaxTransf(const int nsym,
			   const int split,				   
			   int* inp, 
			   int* out,
			   int* ic,
			   int maxtab);

//! Fast Hadamard transform on place 
void FHT(const int n, //!< Log2 of size of transform 
		 int* data	  //!< Input-output data 
		 );

//! Log of a multinomial
double	logmultinom(const int N, const int *a);

//! Range of bit error probabilities
double* range(const int nerr, const int ntot, const double P);

void maxxacc(double & a, const double b);
double maxx(const double a, const double b);
void maxacc(double & a, const double b);


//! Spline interpolation routine
void spline(const int n, const double* x, const double* y,  double* &y2, double yp1=1.e30, double ypn=1.e30);
//! Spline interpolation routine
void splint(const int n ,const double* x, const double* y, const  double* y2, const double xi, double &yi);
//! Linear interpolation routine
void linint(const int n ,const double* x, const double* y, const double xi, double &yi);


int Hweight(int a);

//! Sort a vector of integer
int* Sort(const int n, const int* v);

//! Sort a vector of doubles
int* Sort(const int n, const double* v);

//! Sort a second vector to order the first one
void Reorder(const int n, const int* sort,  int* b);

//! Jacobi algorithm
void jacobi(double *a, int n, double* d, double *v, int *nrot);

//! Cholesky Decomposition of complex Hermitian matrices
double* Cholesky(const double *a,const int n);

//! Singular Value Decomposition of Real matrices
void SVD(double **a, //!< Input complex matrix
		 const int m,  
		 const  int n, 
		 double* w, 
		 double **v
		 );

//! Singular Value Decomposition of Complex matrices (m>=n)
void CSVD(double **a,   //!< Input complex matrix
		  const int m,	//!< Number of rows
		  const int n,	//!< Number of columns
		  double **u,	//!< First nu  of n columns of U 
		  double **v,	//!< First nv  of n columns of V
		  double *s,	//!< Singular values
		  int nu=-1,	//!< Number of computed columns of U
		  int nv=-1,	//!< Number of computed columns of V
		  int p=0		//!< Precomputed inputs
		  );

//! Return the convolution of two real vectors
double* Convolution(const int n1, const double* a1,const int n2, const double* a2, const int rev=0);

//! Modified Bessel functions
double bessi0(double x);
double bessi1(double X);
double bessi(int N, double X);
double bessiN(int N, double X); // Normalized to bessi(0,x)

//! returns log[I0(x)]
double logbessi0(double x);

//! returns log[I1(x)]
double logbessi1(double x);

//! returns atanh(x)
inline double atanh(const double x){return 0.5*log((1.+x)/(1.-x));}



//! Gray mapping
void Gray( const int m, int* map, const int n=1);

//! Collect LLR on constituent symbols and construct LLRs of composite symbol
void Collect(const int tics,	//!< Number of processed composite symbols
			 const int* Input,	//!< LLRs of constituent symbols
			 int* Output,		//!< LLRs of composite symbols
			 const  int ns,		//!< Number of constituent symbols
			 const int* split	//!< Cardinalities of constituent symbol alphabets
			 );

//! Project  LLRs of composite symbol into LLR on constituent symbols 
void Project(const int tics,	//!< Number of processed composite symbols
			   int* inp,		//!< LLRs of composite symbols
			   int* out,		//!< LLRs of constituent symbols
			   const int nsym,	//!< Cardinality of input alphabet
			   const int nlf,	//!< Number of LLRs of composite symbols
			   const int ns,	//!< Number of constituent symbols
			   const int *split,//!< Cardinalities of constituent symbol alphabets				   
			   int maxtab=0,	//!< Number of entries correction table for max*
			   const int *ic=0	//!< Look-up table for max* correction 
			   );

//! Get the i-th bit from a sequence of integers a
inline int getbit(const unsigned int* a,const int i)
{
	return (a[i/32]>>(i%32))&1;
}

//! Set the i-th bit from a sequence of integers a
inline void setbit(unsigned int* a,const int i,const int bit)
{

	if(bit==1)a[i/32] |=  (1<<(i%32));
	else 	  a[i/32] &= ~(1<<(i%32));
}

//! Solve Ax=b with forward substitution
void ForwardSubstitution(const double *A,	//!<Complex lower triangular matrix
								double *b,	//!<Complex known right-hand side
								double *x,	//!<Complex output vector
								int n,		//!<Size of vector (A nxn - b nx1 - x nx1)
								int step	//!<Output values are stored in x every step complex values
								);

/*! Returns determinant of a real matrix */
double Det(const double * a, const int k);

/*! Returns inverse of a generic (real) matrix */
/*! Allocate the buffer to store the matrix */
double * Inverse_gen(const double * A, const int f, double* det=0);

/*! Matrix multiplication of a (real) matrices */
/*! Allocate the buffer to store the resulting matrix */
double * mmult(const double * A, const double * B, const int N, const int K, const int M);

/*! Matrix multiplication of a binary matrices */
/*! Allocate the buffer to store the resulting matrix */
int * mmult(const int * A, const int * B, const int N, const int K, const int M);


/*! Returns inverse of a generic (binary) matrix */
/*! Allocate the buffer to store the matrix */
int * Inverse_bin(const int * A, const int f);


/*! @}*/
void Cholesky_real(double *l, const double *a,const int n);
//! Matrix Inverse for real symmetric
double* Inverse_real(const double *a,const int n);

//! Matrix Inverse for hermitian matrices
double* Inverse(const double *a,const int n);

double LogEuler(const double q);


#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */