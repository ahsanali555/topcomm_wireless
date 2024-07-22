// Output_Manager.h: interface for the Output_Manager class.
//
//////////////////////////////////////////////////////////////////////


#pragma once
#include <stdio.h>

/*! \file  
\brief Declaration of the class Output_Manager
*/
/*! \ingroup SM */
//! Structure that contain the data to be printed with Output_Manager
struct Param
{
	const char* name;	//!< Name of signal used in the file
	const void* ref;	//!< Pointer to Signal	
	const char* format;	//!< format for printing
	int type;			//!< Type of Signal
	int over;			//!< Oversampling of signal
};


/*! \ingroup SM */
//! Manage output of a set of signals to a file. 
/*! The Output Manager samples
the set of signals in a way specified by the user and writes the
results in the corresponding output file.

Each parameter or signal represents a column in the file and a new
line is added to the file for each new sample. The signal is passed
to the object with the method AddSignal(), where the user can
specify the name of the signal, the reference pointer, and optionally
the format to use and an oversampling factor \f$o\f$, which by default is
1.


The signal index \f$i_s\f$ is obtained from current time index \f$t\f$ and
the oversampling factor \f$o_s\f$ through
\f[
i =\lfloor t/o_s\rfloor
\f]

The optional method SetSampling() allows to print an under-sample
version of the signals. The user specifies the sampling period $N$
and initial offset $d$.
\f[
t= d + n\cdot N,
\f]
By default \f$d=0\f$ and \f$N=1\f$, so that all signal samples are printed.


The Run() method, which accepts as single parameter the number of samples \f$n_s\f$ to be printed,
prints \f$n_s\f$ lines to the output file, corresponding to \f$n_s\f$ samples of the associated signals.
The index \f$i_s\f$ of the \f$s\f$-th signal for the \f$n\f$-th sample is then obtained as follows
\f[
    i_s = \frac{d + n\cdot N}{o_s}\;\;\;n=0,\ldots,n_s-1
\f]

 */
class Output_Manager  
{
public:

	//! Default Constructor
	Output_Manager(
		int nmax = 512							//!< Maximum number of parameters
		);	
	
	//! Overload Constructor
	Output_Manager(const char* name,//!< Name of output file
		const bool single = false,				//!< Avoid overwriting of file (auto-renaming)?
		int nmax = 512							//!< Maximum number of parameters
		);	

	virtual ~Output_Manager();

#ifdef WIN32 
	//! Set the name of the output file
	/*! The single flag allow to avoid the overwriting of an existing file. If the flag is on 
	the methods checks if a file exists with the same name. In that case a numbered suffix 
	is added to the provided string to create always a new file.*/
	void SetFile(const char* name, //!< name of file
				const bool single=false //!< Avoid overwriting of file (auto-renaming)?
				);
#endif
	//! Add a signal to the output (integer overload)
	void Add_Signal(const char* name,	//!< Name of signal
		  const int* ref,				//!< Reference
		  const char* format="%d",		//!< Format
		  const int over=1				//!< oversampling of signal
		  );

	//! Add a signal to the output (float overload)
	void Add_Signal(const char* name,	//!< Name of signal
		  const float* ref,				//!< Reference
		  const char* format="%f",		//!< Format for printing
		  const int over=1				//!< oversampling of signal
		  );

	//! Add a signalto the output (double  overload)
	void Add_Signal(const char* name,		//!< Nameof signal
				const double* ref,			//!< Reference
				const char* format="%lf",	//!< Format for printing
				const int over=1			//!< Oversampling of signal
				);

	//! Add a complex signal to the output
	void Add_ComplexSignal(const char* name,	//!< Name of signal
				const double* ref,				//!< Reference
				const char* format="%lf",		//!< Format for printing
				const int over=1				//!< Oversampling of signal		
				);

	//! Print in the first column the time axis
	void Add_Time(const double dt=1. //!< Add in the first column the time axis
	);

    void Header(FILE*file = NULL);
	void Run(const int ntics);
	void RunSingle(FILE* file = NULL);
	void Display(FILE* file = stdout) const;
	bool Check(const int ntics);

	//! Set the sampling period and offset of signals (default to one)
	void SetSampling(const int offin,const int nsampin);


private:
	int Add_Signal(const char* name,	// name
				  const char* format, // format
				  const void* ref,	// reference
				  const int type,	// reference
				  const int over	// reference
				  );
	int nmax;
	int nparam;
	Param **param;		// List of Signals
	bool time;	        // time axis is reported
	FILE* stream;		// stream pointer
	bool first;
	double dt;	// Delta associated to one tics of time
	int nclock;
	int off;
	int nsamp;
};

