#pragma once
#include "Demodulator.h"
/*! \file
\brief Declaration of class Tuning_Centroids */
/*! \ingroup Modems
\brief Centroid tuning algorithm.
Estimates the position of centroids (average positions of observed constellation points) and their variance. 

The estimate is performed by averaging the provided observation samples 
conditioned on the knowledge of transmitted constellation labels, which are also passed to this
block.

\author Guido Montorsi
*/
class Tuning_Centroids
{
public:
	Tuning_Centroids();
	~Tuning_Centroids(void);

	//! Set the main parameters of class
	void SetParameters(
		const Modulator*,     //!< Pointer to reference modulator, specifying the reference constellation			
		const double alpha=0. //!< Bandwidth of estimation one pole filter (default to zero bandwith (accumulator))
		);

	//! Reset the centroid statistic
	void Reset();

	//! Update the statistic based on the sequence of transmitted indexes
	/*! The sequence of transmitted indexes may  also be an estimate. A "-1" in the sequence
	of transmitted indexes indicates that the index is not available and the corresponding observed 
	sample is dropped from the statistic. */
	void Update_Positions(
		const int nsymb,	//!< [in] Number of observed samples
		const double* modRX,//!< [in] Observed samples 
		const int* ssymb	//!< [in] Indexes of transmitted points (before mapping)
		);

	//! Update a Modulator and a Demodulator with  the current centroid statistic
	/*! The constellation set used by the Modulator is updated with the current estimated centroid position.
	The variance used by the Demodulator for the LLR computation is updated with the current estimated
	variance of centroids.
	*/
	void Set_Demodulator(
		Modulator* refmod, //!< Modulator to be updated.
		Demodulator* Demod  //!< Demodulator to be updated
		) const;

	//! Return the current centroid estimation
	double Get_Centroids(
		double* mod,   //!< Centroid positions [output, vector of M values]
		double* var = 0, //!< Variances of each centroid [optional vector of M values]
		double*ISI = 0   //!< Average variance of centroids [optional pointer to scalar]
		) const;

	//! Display current centroid statistic
	void Display_Centroids(
		FILE * file=stdout,		//!< Output file
		const bool verbose=false //!< Flag for verbose output
		) const;

	//! Get the currently estimated SNR [dB]
	/*! The ratio between the average energy of the centroids and the average variance */
	double GetSNR() const;

private:
	const Modulator* mod;
	double alpha;
	double* centroids;
	double* vqm;
	__int64* stat;
	int M;

};
