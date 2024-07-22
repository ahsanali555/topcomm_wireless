#pragma once
#include <math.h>
#include <stdio.h>
class Frame_Sync
{
public:

	Frame_Sync(void);

	//! Set the reference TX system and the number of SCCC decoder iterations 
	void SetParameters(
		const int FL,
		const int H,
		const double *Header,
		const int ncand = 10,
		const int coh = 16,
		const double afoff = 0.01);

	void Reset();

	//! Run frame sync algorithm and frequency offset estimation
	/*! Returns the number of delivered frames. */
	int Run(
		const int ntics,
		const double*inp,
		double*out=0 //!< If output buffer is provided, aligned frame is provided at the output
		);



	~Frame_Sync(void);
	void Display(FILE * file=stdout) const;

	/**\{ */
	int sync;
	double foff;   // Current frequency estimation
	int pointer;   // Pointer to next delivered frame
	int ini;       // Status of sync

private:
	int FL;
	int H;
	double phcur;    // Current phase
	double* Header;	 // Header 
	double* best;    // Candidates
	double* line;    // shift register 
	double* buffer;  // To store aligned frame

	int coh;		
	double afoff;
	int ncand;
	int time ;
	int c;
	int jmax;       // Best candidate in the list
	int jmin;		// Worst candidate in the list
};
