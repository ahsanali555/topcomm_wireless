// Parameter_Manager.h: interface for the Parameter_Manager class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdio.h>
/*! \file 
\brief Declaration of the class Parameter_Manager
*/
/*! \ingroup SM 
\brief Manages a set of simulation parameters.

This utility class is used in the main programs to manage the input of parameters from a
user file and to update them to produce automatically multiple simulations with variable
system parameters.

The parameters are read from a file and the possible values of 
the parameter can be specified as 
- constants
- set of possible values with starting, steps and maximum value
- Set of values
Each parameter have a name string , that is reported in the line preceding the
specification of its (possible) values.

To specify a set of uniformly spaced possible values, use the convention 
<min>:<max>:<step>
This will generate a set of simulations where the parameter takes the all the 
values between <min> and <max> incremeting by <step>.


To specify a generic set of possible values, use the convention 
<a>|<b>|<c>|...
This will generate a set of simulations where the parameter takes the value <a>, then <b>, 
and so on.

If more parameters are specified through sets, ALL the possible simulation
configurations will be generated, corresponding to the cartesian product of the specified
sets.
The order used to change the parameters depends on the order the parameters are read by
the main file.

See the manual for more details and examples.

*/

class Parameter_Manager  
{
public:
	Parameter_Manager(int np=40);
	virtual ~Parameter_Manager();


	//! Add a new parameter and get its values from a file (overload for integer par.) 
	bool Add_Parameter_From_File(
		FILE* file,		//!< Pointer to the file containing the parameter values.
		char* name,		//!< Name of the string to be searched and of the parameter.
		int* var,		//!< Pointer to the variable containing the parameter(s)
		char* = "%d"	//!< Format for reading the values from the file
		);

	//! Add a new parameter and get its values from a file (overload for floats par.) 
	bool Add_Parameter_From_File(
		FILE* file,		//!< Pointer to the file containing the parameter values.
		char* name, 	//!< Name of the string to be searched.
		float* var,		//!< Pointer to the variable containing the parameter(s)
		char* = "%f"	//!< Format for reading the values from the file
		);

	//! Add a new parameter and get its values from a file (overload for double par.)  
	bool Add_Parameter_From_File(
		FILE* file,		//!< Pointer to the file containing the parameter values.
		char* name,		//!< Name of the string to be searched.
		double* var,	//!< Pointer to the variable containing the parameter(s)
		char* = "%lf"	//!< Format for reading the values from the file
		);

	//!Add a new parameter and get its values from a file (overload for string par.)
	bool Add_Parameter_From_File(
		FILE* file,			//!< Pointer to the file containing the parameter values.
		char* name,			//!< Name of the string to be searched.
		char* var,			//!< Pointer to the variable containing the parameter(s)
		char* = "%[^|^\n^,]"	//!< Format for reading the values from the file
		);

	//! Add a new parameter to the list of managed parameters.
	void Add_Parameter(
		char* name, //!< Name of the string to be searched and of the parameter.
		void* var, 	//!< Pointer to the variable containing the parameter(s)
		void* val, 	//!< Pointer to the vector values assumed by the parameter
		int type,	//!< Type of parameter (0: integer, 1: float, 2: double, 3: strings).
		int size=1  //!< Size of value vector
		);

	//!< Set the flags to control where the parameter is shown (1: List, 2: online, 3: filename)
	void SetOutflags(const int out = 0xff);

	//!< Get a file name
	char* GetFileName() const;


	//! Add a new constant parameter (overload for integers) 
	void Add_Constant(
		char* name, 
		int* var
		);

	//! Add a new constant parameter  (overload for floats) 
	void Add_Constant(
		char* name, 
		float* var
		);

	//! Add a new constant parameter  (overload for doubles) 
	void Add_Constant(
		char* name, 
		double* var
		);

	//! Add a new constant parameter (overload for strings) 
	void Add_Constant(
		char* name, 
		char* var
		);
	
	//! Update the parameter list. Return false if no updating is possible 
	bool Update_Parameters();

	//! Print the parameters list 
	void Print_Parameters(FILE* =stdout) const;

	//! Print the parameters on a line spaced by tabs
	/*! When the header flag is set to true, the list of parameter 
	names is printed instead of the list of 
	parameters values. This is useful at the beginning of simulation for the header.*/
	void Print_Parameters_on_Line(
		FILE *file=stdout,	//!< Output stream (default to stdout)
		bool header=false	//!< Print the names or the parameter values?
		) const;

	//! Print the list of the changed parameters in the last simulation run
	void Print_Changed(FILE *file=stdout) const;

	bool OuterLoop() const;

	bool error;
	int changedpar;		// Changed parameter

private:
	void Add_Constant(char* name, void* var, int type);
	bool Add_Parameter_From_File(FILE* file,char* name, void* var,int type, char*);

	void	**parameters;	// array of Pointer parameters
	int		*types;			// Array of types
	void	**values;		// array of pointer to lists
	char	**names;		// Array of names
	int		*sizes;			// Size of list of values
	int		*current;		// current points
	int	    *outflag;		    // flag for printing it
	int		npar;				// Current number of parameters
	bool find_line(FILE* file,char* name);
	int maxstring;
};
