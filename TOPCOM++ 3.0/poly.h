#pragma once
/** \cond INTERNAL */
/*!
\ingroup Poly
 \file 
 \brief Functions to manage polynomials with complex coefficients and rational functions.

  This set of functions are used by the interface functions of the Filter.

  \see 
  		- Chebyshev()
		- Elliptic()
		- Bessel()
		- Butterworth()

*/
/*@{*/

  
//! Derivative a polynomial.
double* DerPoly(const int n, const double* pol);

//! Indefinite integral a polynomial.
double* IntPoly(const int n,  const double* pol);

//! Product of two polynomials.
double* ProdPoly(const int n1, const double* pol1,const int n2, const double* pol2);

//! Sum of two polynomials.
double* SumPoly(const int n1, 
		   const double* pol1,
		   const int n2, 
		   const double* pol2);

//! Square absolute value of a polynomial.
double* AbsPoly(const int n1, const double* pol1);

//! Division of two polynomial.
double* DivPoly(const int n1, const double* pol1,const int n2, const double* pol2);

//! Decomposition of rational function into simple ratios.
double* DecompRatio(const int n, const double* num,const int d, const double* den, const int *mul=0);

//! Derivative a a rational function.
void	DerRatio(const int n, 
				 const double* num, 
				 const int d,
				 const double* den,
				 int	 &nd, 
				 double* &numd, 
				 int	 &dd,
				 double* &dend);

//! Bilinear transformation of a rational function.
void	Bilinear(const int n, 
				 const double* num, 
				 const int d,
				 const double* den,
				 int	 &nd, 
				 double* &numd, 
				 int	 &dd,
				 double* &dend);

//! Evaluate a polynomial in a point 
void EvalPoly(const int n, 
			  const double* pol, 
			  const double xR,
			  const double xI,
			  double &yR,
			  double &yI
			  );


//! Sum two rational functions
void SumRatio(
		const int n1, 
	    const double* num1, 
	    const int d1,
	    const double* den1,
		const int n2, 
	    const double* num2, 
	    const int d2,
	    const double* den2,
		int		&nd, 
	    double* &numd, 
	    int		&dd,
	    double* &dend
		);

/*#include "Filter.h"
@}*/