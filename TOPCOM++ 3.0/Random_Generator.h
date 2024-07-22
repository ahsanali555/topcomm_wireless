// Random_Generator.h: interface for the Random_Generator class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#ifdef CTOPCOM
#include "ctopcom.h"
#endif
/*! \file
\brief Declaration of classes Random_Generator, Gaussian_Generator, Rayleigh_Generator, Rice_Generator, LogNormal_Generator.
*/

/*! \ingroup Random_Sources 
\brief  Generate arbitrary random numbers.

The Random_Generator class generates arbitrary random variables using the rejection method.
User can specify the desired distribution through 
	- a table 
	- a pointer to the function describing the distribution,
	- a pointer to the function returning the inverse cumulative function

\author Guido Montorsi
*/
class Random_Generator  
{
public:
	//! Set the seed used to inizialize the random generator
	Random_Generator(int seed=1		//<! Seed used to inizialize the random generator (default:1)
	);
	virtual ~Random_Generator();	

	//! Set the distribution and its support
	void SetDistribution(
			double f(double),    //!< Distribution function (cdf).
		   double xmin=-100.,	//!< Minimum abscissa of the distribution.
		   double xmax= 100.,	//!< Maximum abscissa of the distribution.
		   int npoints=1000		//!< Number of points of the cdf.
		   );
	
	//! Set the inverse cumulative function
	void SetInvCumulative(
		double f(double)	//!< A function that return the inverse cumulative function of the desired distribution
		); 
	
	//! Set the seed of the generation
	void SetParameters(const int seed=1		//!< Seed for the generation of uniform random variable (default:1).
			);
	
	//! Set the distribution type to uniform
	void SetUniform(const bool u=true){uniform=u;};
	
	void SetTikhonov(const double alpha, const int nsamp=1000);


	//! Generate the random variables
	int Run(const int tics,			//!< Number of output random variables
			double* out				//!< Output r.v. vector
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,			//!< Number of output random variables
			cmplx* out				//!< Complex Output r.v. vector
		)
	{
	return(Run(tics,(double*) out));
	}
#endif


private:
	int seed;			//!< Seed of random generator
	bool uniform;		//!< If unform=true, a uniformly distributed variable in (0,1) is generated
	double* invcum;		//!< Inverse cumulative
	int npoints;		//!< Number of points
	double (*f)(double);	//!< Cumulative distribution function
	double alpha;		//!< Used for Tik. distribution
};


/*!\ingroup Random_Sources 
\brief  Generate random numbers with Gaussian distribution.

  \f[
f(y|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}{\rm e}^{-\frac{(y-\mu)^2}{2\sigma^2}}
\f]
The user specifies the random generator 
seed and the parameters (\f$\sigma\f$ and \f$\mu\f$) of the distribution.

For an example of its use see e.g. the test program "test_random_generator.cpp".

\author Gabriella Bosco
*/
class Gaussian_Generator  
{
public:
	//! Set the seed used to inizialize the random generator
	Gaussian_Generator(int seed=1		//<! Seed
	);
	virtual ~Gaussian_Generator();	

	//! Set the parameters of the Gaussian distribution
	void SetParameters(const int seed=1,		//<! Seed
					const double sigma=1.,		//<! Parameter \f$\sigma\f$
					const double mu=0.			//<! Parameter \f$\sigma\f$
					);

	//! Set the seed of the generator
	void SetSeed(const int seed=1		//<! Seed
		);	

	//! Generate the random variables
	int Run(const int tics,			//!< Number of output random variables
			double* out				//!< Output r.v. vector
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,			//!< Number of output random variables
			cmplx* out				//!< Complex Output r.v. vector
		)
	{
	return(Run(tics,(double*) out));
	}
#endif

private:
	double sigma,	//!< Parameter \f$\sigma\f$ of the distribution 
		mu;			//!< Parameter \f$\mu\f$ of the distribution 
	double x2;			//!< Temporary variable
	bool isready;	//!< Is the r.v. ready?
	int seed;			//!< Seed of random generator
};

/*!\ingroup Random_Sources 
\brief  Generate random numbers with Log Normal distribution.

\f[
f(y|\mu,\sigma)=\frac{1}{x\sigma\sqrt{2\pi}}{\rm e}^{-\frac{(\log{(x)}-\mu)^2}{2\sigma^2}}.
\f]
with  \f$\mu\f$ and \f$\sigma\f$ defined by the user through the method SetParameters().
A random variable \f$y\f$ with lognormal distribution with parameters \f$\mu\f$ and \f$\sigma\f$
 is generated from a
random variable \f$x\f$ with Gaussian distribution, \f$N(0,1)\f$, using the following transformation:
\f[
y=\exp\{\sigma x+\mu\}.
\f]
A tic of the Run() method generates one random number with Lognormal distribution.

For an example of its use see e.g. the test program "test_random_generators.cpp".

\author Gabriella Bosco
*/


class LogNormal_Generator  
{
public:
	//! Set the seed used to inizialize the random generator
	LogNormal_Generator(int seed=1		//<! Seed 
	);
	virtual ~LogNormal_Generator();	
	//! Set the parameters of the LogNormal distribution
	void SetParameters(
		const int seed=1,		//<! Seed 
		const double sigma=1.,		//<! Parameter \f$\sigma\f$ 
					const double mu=0.			//<! Parameter \f$\mu\f$ 
					);	

	//! Set the seed of the generator
	void SetSeed(const int seed=1		//<! Seed 
		);		
	

	//! Generate the random variables
	int Run(const int tics,			//!< Number of output random variables
			double* out				//!< Output r.v. vector
		);	

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,			//!< Number of output random variables
			cmplx* out				//!< Complex Output r.v. vector
		)
	{
	return(Run(tics,(double*) out));
	}
#endif

private:
	double sigma,	//!< Parameter \f$\sigma\f$ of the distribution 
		mu;			//!< Parameter \f$\mu\f$ of the distribution 
	Gaussian_Generator* GG;	//!< Gaussian r.v. generator
};

/*!\ingroup Random_Sources 
\brief  Generate random numbers with Rice distribution.

\f[
f(y|\mu,\sigma)=\frac{y}{\sigma^2}{\rm
e}^{\frac{-(y^2+\mu^2)}{2\sigma^2}}I_0\left(\frac{y\mu}{\sigma^2}\right)\;\;\;y>0
\f]
with  \f$\mu>0\f$ and \f$\sigma\f$ defined by the user through the method SetParameters().
\f$I_0\f$ is a modified Bessel function of the first kind.

A random variable \f$y\f$ with Rice distribution with parameters \f$\mu\f$ and \f$\sigma\f$ 
is generated from two
random variables \f$x_1\f$ and \f$x_2\f$ with Gaussian distribution \f$N(0,1)\f$
 using the following transformation:
\f[
y=\sqrt{(\sigma x_1+\mu)^2+(\sigma x_2)^2}.
\f]
A tic of the Run() method generates one random number with Rice distribution.

For an example of its use see e.g. the test program "test_random_generator.cpp".

\author Gabriella Bosco
*/
class Rice_Generator  
{
public:
	//! Set the seed used to inizialize the random generator
	Rice_Generator(int seed=1		//<! Seed 
	);
	virtual ~Rice_Generator();	
	
	//! Set the parameters of the Rice distribution
	void SetParameters(const int seed=1,		//<! Seed (
					const double sigma=1.,		//<! Parameter \f$\sigma\f$ 
					const double mu=0.			//<! Parameter \f$\mu\f$
					);	
	//! Set the seed of the generator
	void SetSeed(const int seed=1		//<! Seed 
		);		
	
	//! Generate the random variables
	int Run(const int tics,			//!< Number of output random variables
			double* out				//!< Output r.v. vector
		);	

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,			//!< Number of output random variables
			cmplx* out				//!< Complex Output r.v. vector
		)
	{
	return(Run(tics,(double*) out));
	}
#endif

private:
	double sigma,	//!< Parameter \f$\sigma\f$ of the distribution 
		mu;			//!< Parameter \f$\mu\f$ of the distribution 
	Gaussian_Generator* GG;	//!< Gaussian r.v. generator
};

/*!\ingroup Random_Sources 
\brief  Generate random numbers with Rayleigh distribution.

\f[
f(y|\sigma)=\frac{y}{\sigma^2}{\rm e}^{-\frac{y^2}{2\sigma^2}}\;\;\;y>0
\f]
with  \f$\sigma\f$ defined by the user through the method SetParameters().
The random variable f$y\f$ with Rayleigh distribution with parameter \f$\sigma\f$ is generated from a
random variable \f$u\f$ with uniform distribution  in (0,1) using the following transformation:
\f[
y=\sigma\sqrt{-2\log{(u)}}.
\f]
A tic of the \Run() method generates one random number with Rayleigh distribution.

\author Gabriella Bosco
*/
class Rayleigh_Generator  
{
public:
	//! Set the seed used to inizialize the random generator
	Rayleigh_Generator(int seed=1		//<! Seed 
	);
	virtual ~Rayleigh_Generator();	
	//! Set the parameters of the Rayleigh distribution
	void SetParameters(const int seed=1,		//<! Seed
					const double sigma=1.		//<! Parameter \f$\sigma\f$ 
					);	
	//! Set the seed of the generator
	void SetSeed(const int seed=1		//<! Seed 
		);		
	
	//! Generate the random variables
	int Run(const int tics,			//!< Number of output random variables
			double* out				//!< Output r.v. vector
		);

#ifdef CTOPCOM
	//! Overload of Run with complex signals
	int Run(const int tics,			//!< Number of output random variables
			cmplx* out				//!< Complex Output r.v. vector
		)
	{
	return(Run(tics,(double*) out));
	}
#endif

private:
	double sigma;	//!< Parameter \f$\sigma\f$ of the distribution 
	int seed;			//!< Seed of random generator
};


