// SCCC_Decoder.h: interface for the SCCC_Decoder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "SCCC_Encoder.h"
#include "SISO_Decoder.h"
#include "Puncturer.h"

/*! \file
\brief Declaration of class SCCC_Decoder.
 */
/*! \ingroup Codecs
\brief Binary Serial Concatenated Convolutional Code Decoder.

The SCCC_Decoder class is associated to
the correspondent SCCC_Encoder. 

Through the method SetParameters the user specifies
- A pointer to the correspondent SCCC\_Encoder that describes the code
- The number of iterations \f$n_I\f$ of the decoder
- The size of updating and training window \f$T_W\f$ and \f$U_W\f$ of the embedded SISO decoders. 
(Optional, default to automatic computation)
- The scaling factor \f$f\f$ of the input quantized LLR. (Optional, default to 8, which gives ideal performances)

Optionally, with the method SetStop the user can activate or deactivate a stopping rule of the decoder. 

The Run method takes as input the quantized LLR of the coded bits (,
after a possible depuncturing, which inserts zero values for those bit that were not transmitted, 
the iterative decoder is performed for \f$n_I\f$ iterations. 
The output buffer is then filled with the last computed  LLRs of the information
bits.


	For an example of its use see e.g. the test program "test_SCCC.cpp".

\author Guido Montorsi
*/
class SCCC_Decoder  
{
public:
	SCCC_Decoder();
	virtual ~SCCC_Decoder();
	void SetParameters(const SCCC_Encoder* refcodin,	//!< Reference SCCC encoder
						const int niter,	//!< Number of iterations
						const int nini=-1,		//!< Training window for SISOs (-1=automatic computation)
						const int ngroup=-1,	//!< Updating window for SISOs (-1=automatic computation)
						const double factor=8.	//!< factor for maxx operation
				);
	void TuneParameters(const int niter);

	//! Activate or deactivate a stopping rule for the decoder
	/*! The stopping rule is based on the comparison of the LLRs signs at the input 
	and at the output of the outer SISO processor. If the signs agree for all bits the 
	iterative decoding is stopped. The stopping rule is activated by default.*/
	void SetStop(const bool stopin=true){stop=stopin;}

	//! Activate or deactivate the external reset of extrinsic information
	/*! Usually extrinsic information is reset to zero at the beginning of
	the first iteration. 	With this method this behavior can be changed, 
	so that the extrinsic information is
	reset only with explicit calls to the ResetExtrinsic method. */
	void SetExternReset(const bool extres=true){this->extres=extres;}

	//! Reset the extrinsic information to zero
	void ResetExtrinsic();

	//! Run the decoder
	int Run(const int tics,					//!< Number of decoded blocks
		const int* inp,						//!< Input LLR on coded bits
		int* out,							//!< Output LLR on inf  bits
		const int* inpin =0,				//!< Input  LLR on inf bits [Optional]
		int* codout=0,						//!< Output LLR on coded bits [Optional]
		const int* data=0					//!< Reference data for genie aided stop [Optional]
		);

	// Display the main configuration parameters of the block
	void Display(FILE* f=stdout);

//	int ave;
//	int* ave2;

	//! Flag for activate stopping rule
	bool stop;
	//! Flag to indicate if resetting of extrinsic is performed at the beginning of operations
	bool extres;			

	int *ext1;
	int *ext2;

	//! Outer SISO decoder
	SISO_Decoder * oSISO;
	//! Inner SISO decoder 
	SISO_Decoder * iSISO;

	friend class  RX_CCSDS_SCCC;
private:
	int niter;
	const SCCC_Encoder *refcod;		//!< Reference to the correspondent encoder.
	Puncturer *puni,*puno;
	int* buffllr1;
	int* buffllr2;
	int* buffllr3;
	bool resetext;
	// For statistic on LLRs
};
