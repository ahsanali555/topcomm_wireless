// Modulator.h: interface for the Modulator class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <stdio.h>
#ifdef CTOPCOM
#include "ctopcom.h"
#endif

/*! \file 
\brief Declaration of class Modulator and its interfaces
*/
/*! \ingroup Modems
\brief Generic memoryless modulator.
\author Guido Montorsi

The Modulator class implements the modulation function for any kind of one or two dimensional 
digital modulation scheme.
The Modulator maps each symbol (composed by \f$m\f$ bits) onto a set of 
\f$M=2^m\f$ constellation points, that are transmitted sequentially over the channel. 
Each signal \f$s_k(t)\f$ is represented by a point \f${\bf s}_k\f$ in the one or 
two dimensional Euclidean space.
The user can specify:
- the constellation points
- the mapping between binary input and one or two dimensional constellation points.

\see Several interface functions are provided to generate standard modulations:
- PSK_Modulator
- PAM_Modulator
- QAM_Modulator
- MultiPSK_Modulator
- APSK_Modulator
- Gaussian_Modulator
*/

class Modulator  
{
public:
Modulator();
virtual ~Modulator();
//! Set the main parameters of the modualtor
int SetParameters(
	const int m,			//!< Constellation cardinality (log2)
	const double* constel,	//!< Constellation
	const int*    mapping=0,//!< Mapping of bits to constellation points (default to identity mapping)
	const bool    one=false	//!< Flag for one dimensional constellation
	);

//! Run  the modulator
void Run(const int ntics,	//!< Number of input symbols
		const int* inp,		//!< Input bits (or symbols)
		double* out,		//!< Output symbols
		const bool bin=true,	//!< If bin=true, the bits (instead of the symbols) are input to the modulator
		int* symbols=0	//! optional buffer to store signal index sequence
		);
//! Returns the PAPR ratio between peak and average power of constellation
double GetPAPR(
	const bool db=true	//!< returns the value in dB?
	) const;


#ifdef CTOPCOM
//! Overload of Run with complex signals
//! Run  the modulator
void Run(const int ntics,	//!< Number of input symbols
		const int* inp,		//!< Input bits (or symbols)
		cmplx* out,			//!< Complex Output symbols
		const bool bin=true,	//!< If bin=true, the bits (instead of the symbols) are input to the modulator
		int* symbols=0	//! optional buffer to store signal index sequence
		)
{
Run(ntics,inp,(double*) out,bin,symbols);
return;
}
	
#endif

// Friends functions
friend
Modulator* PSK_Modulator(const int m,
        const int* mapping,
        const double avenergy
        );

friend
Modulator* PAM_Modulator(const int m,
        const int* mapping,
        const double avenergy
        );

friend
Modulator* QAM_Modulator(const int m,
        const int* mapping,
        const double avenergy
        );


friend
Modulator* MultiPSK_Modulator(const int Ncircles,
        const int* n,
        const double* rho,
        const double* phi,
        const int* mapping
        );

friend
Modulator* DVBS2_APSK_Modulator(const int m,
        const double spec_eff
        );


friend
Modulator* DVBSX_Modulator(
        const int index
        );



friend
Modulator* Nonuniform_QAM_Modulator(const int m, const int a);


friend
Modulator* xDSL_Modulator(const int m);



friend
Modulator* Gaussian_Modulator(
        const int m,
        const int* mapping
        );


friend
Modulator* Read_from_File(const char* namefile, const int skip=0);

Modulator * Get_Constellation(const int m, char * fileconst, const int index, const int oneD);

friend
Modulator *Get_Constellation(
	const int m,        //!< Modulation efficiency (bits/symbol)
	char * fileconst,   //!< Name of file or ("CCSDS","QAM", DVBS2")
	const int index=0,   //!< Additional index to specify modulation set
	const int oneD=0	//!< Additional index to specify if One D constellation is loaded
	);
void Display(FILE* = stdout);
void Normalize(const double aven=1, const bool peak=false);

friend class Demodulator;
friend class MIMO_Demodulator;
friend class TCM_Encoder;
double* constellation;  //!< Points of constellation
int m;					//!< \f$\log_2(M)\f$
bool	one;			//!< Flag for mondimensional constellation
int M;					//!< Constellation cardinality
int*    mapping;		//!< Mapping of bits to constellation point
private:

};

/* Gray Mapping (up to 1024)*/
static int gray[1024]=
{0,1,3,2,7,6,4,5,15,14,12,13,8,9,11,10,
31,30,28,29,24,25,27,26,16,17,19,18,23,22,20,21,
63,62,60,61,56,57,59,58,48,49,51,50,55,54,52,53,
32,33,35,34,39,38,36,37,47,46,44,45,40,41,43,42,
127,126,124,125,120,121,123,122,112,113,115,114,119,118,116,117,
96,97,99,98,103,102,100,101,111,110,108,109,104,105,107,106,
64,65,67,66,71,70,68,69,79,78,76,77,72,73,75,74,95,94,92,93,88,89,91,90,80,81,83,82,87,86,84,85,
255,254,252,253,248,249,251,250,240,241,243,242,247,246,244,245,224,225,227,226,231,230,228,229,
239,238,236,237,232,233,235,234,192,193,195,194,199,198,196,197,207,206,204,205,200,201,203,202,
223,222,220,221,216,217,219,218,208,209,211,210,215,214,212,213,128,129,131,130,135,134,132,133,
143,142,140,141,136,137,139,138,159,158,156,157,152,153,155,
154,144,145,147,146,151,150,148,149,191,190,188,
189,184,185,187,186,176,177,179,178,183,182,180,181,160,161,163,162,167,166,164,165,175,174,172,173,168,169,171,
170,511,510,508,509,504,505,507,506,496,497,499,498,503,502,500,501,480,481,483,482,487,486,484,485,495,494,492,
493,488,489,491,490,448,449,
451,450,455,454,452,453,463,462,460,461,456,457,459,458,479,478,476,477,472,473,475,474,464,465,467,466,471,470,
468,469,384,385,387,386,391,390,388,389,399,398,396,397,392,393,395,394,415,414,412,413,408,409,411,410,400,401,
403,402,407,406,404,405,447,446,444,445,440,441,443,442,432,433,435,434,439,
438,436,437,416,417,419,418,423,422,420,421,431,430,428,429,424,425,427,426,256,257,259,258,263,262,260,261,271,
270,268,269,264,265,267,266,287,286,284,285,280,281,283,282,272,273,275,274,279,278,276,277,319,318,316,317,312,
313,315,314,304,305,307,306,311,310,308,309,288,289,291,290,295,294,292,293,
303,302,300,301,296,297,299,298,383,382,380,381,376,377,379,378,368,369,371,370,375,374,372,373,352,353,355,354,
359,358,356,357,367,366,364,365,360,361,363,362,320,321,323,322,327,326,324,325,335,334,332,333,328,329,331,330,
351,350,348,349,344,345,347,346,336,337,339,338,343,342,340,341,1023,1022,1020,
1021,1016,1017,1019,1018,1008,1009,1011,1010,1015,1014,1012,1013,992,993,995,994,999,998,996,997,1007,1006,1004,
1005,1000,1001,1003,1002,960,961,963,962,967,966,964,965,975,974,972,973,968,969,971,970,991,990,988,989,984,985,
987,986,976,977,979,978,983,982,980,981,896,897,899,898,903,902,900,901,
911,910,908,909,904,905,907,906,927,926,924,925,920,921,923,922,912,913,915,914,919,918,916,917,959,958,956,957,
952,953,955,954,944,945,947,946,951,950,948,949,928,929,931,930,935,934,932,933,943,942,940,941,936,937,939,938,
768,769,771,770,775,774,772,773,783,782,780,781,776,777,779,778,799,798,796,
797,792,793,795,794,784,785,787,786,791,790,788,789,831,830,828,829,824,825,827,826,816,817,819,818,823,822,820,
821,800,801,803,802,807,806,804,805,815,814,812,813,808,809,811,810,895,894,892,893,888,889,891,890,880,881,883,
882,887,886,884,885,864,865,867,866,871,870,868,869,879,878,876,877,872,873,
875,874,832,833,835,834,839,838,836,837,847,846,844,845,840,841,843,842,863,862,860,861,856,857,859,858,848,849,
851,850,855,854,852,853,512,513,515,514,519,518,516,517,527,526,524,525,520,521,523,522,543,542,540,541,536,537,
539,538,528,529,531,530,535,534,532,533,575,574,572,573,568,569,571,570,560,
561,563,562,567,566,564,565,544,545,547,546,551,550,548,549,559,558,556,557,552,553,555,554,639,638,636,637,632,
633,635,634,624,625,627,626,631,630,628,629,608,609,611,610,615,614,612,613,623,622,620,621,616,617,619,618,576,
577,579,578,583,582,580,581,591,590,588,589,584,585,587,586,607,606,604,605,
600,601,603,602,592,593,595,594,599,598,596,597,767,766,764,765,760,761,763,762,752,753,755,754,759,758,756,757,
736,737,739,738,743,742,740,741,751,750,748,749,744,745,747,746,704,705,707,706,711,710,708,709,719,718,716,717,
712,713,715,714,735,734,732,733,728,729,731,730,720,721,723,722,727,726,724,
725,640,641,643,642,647,646,644,645,655,654,652,653,648,649,651,650,671,670,668,669,664,665,667,666,656,657,659,
658,663,662,660,661,703,702,700,701,696,697,699,698,688,689,691,690,695,694,692,693,672,673,675,674,679,678,676,
677,687,686,684,685,680,681,683,682};

//{
//	0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13, 8, 9, 11, 10, 
//	31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21,
//	63, 62, 60, 61, 56, 57, 59, 58, 48, 49, 51, 50, 55, 54, 52, 53, 
//	32, 33, 35, 34, 39, 38, 36, 37, 47, 46, 44, 45, 40, 41,	43, 42, 
//	127, 126, 124, 125, 120, 121, 123, 122, 112, 113, 115, 114, 119, 118, 116, 117,
//	96, 97, 99, 98, 103, 102, 100,	101, 111, 110, 108, 109, 104, 105, 107, 106, 
//	64, 65, 67, 66, 71, 70, 68, 69, 79, 78, 76, 77, 72, 73, 75, 74, 
//	95, 94, 92, 93, 88, 89, 91, 90, 80, 81, 83, 82, 87, 86, 84, 85
//};

static int Mapping8PSK[8]={0,1,3,2,7,6,4,5}; // MSB protected	
static int Mapping8PSKinv[8]={0,7,3,4,1,6,2,5}; // LSB protected


/* 128 Mapping constructed by Ten Brink */
static int Mapping128QAMGray1[128]=
{
18,30,54,42,102,90,66,78,19,31,55,43,103,91,67,79,21,33,57,45,105,93,
69,81,20,32,56,44,104,92,68,80,25,37,61,49,109,97,73,85,24,36,60,48,
108,96,72,84,22,34,58,46,106,94,70,82,23,35,59,47,107,95,71,83,17,29,
53,41,101,89,65,77,16,28,52,40,100,88,64,76,11,3,8,0,115,123,112,120,
10,2,9,1,114,122,113,121,26,38,62,50,110,98,74,86,27,39,63,51,111,99,
75,87,12,4,15,7,116,124,119,127,13,5,14,6,117,125,118,126
};

static int Mapping32QAMGray[32]=
{
12,0,6,1,11,10,5,4,18,22,24,29,17,16,23,28,
13,9,7,2,14,15,8,3,19,31,25,30,20,21,26,27
};


/*! \ingroup Interface
\brief Generate a PSK modulator
*/

Modulator* PSK_Modulator(const int m, //!< Number of modulation bits
    const int* mapping = 0, //!< Signal mapping (0: Gray mapping, 1 is natural mapping)
    const double avenergy = 1. //!< Average signal energy 
    );
/*! \ingroup Interface
\brief Generate a PAM modulator
*/

Modulator* PAM_Modulator(const int m, //!< Number of modulation bits
    const int* mapping = 0, //!< Signal mapping (0: Gray mapping, 1 is natural mapping)
    const double avenergy = 1. //!< Average signal energy 
    );
/*! \ingroup Interface
\brief Generate a QAM modulator
*/

Modulator* QAM_Modulator(const int m, //!< Number of modulation bits
    const int* mapping = 0, //!< Signal mapping (0: Gray mapping, 1 is natural mapping)
    const double avenergy = 1. //!< Average signal energy
    );

/*! \ingroup Interface
\brief Generate a circular PSK modulator

This interface is used to generate a multi-circular PSK constellation, 
composed by $N$ concentric circular PSK constellations.
*/

Modulator* MultiPSK_Modulator(const int Ncircles, //!< Number of PSK circles 
    const int* n, //!< Number of points in each circle 
    const double* rho, //!< Radius of each circle 
    const double* phi, //!< Phase offset of each circle
    const int* mapping = 0 //!< Signal mapping (default: Natural mapping)
    );
/*! \ingroup Interface
\brief Generate an APSK modulator as in the DVB-S2 standard
*/

Modulator* DVBS2_APSK_Modulator(const int m, //!< Modulation efficiency (available values: 4,5 or 6)
    const double spec_eff = -1 //!< Spectral efficiency
    );


/*! \ingroup Interface
\brief Generate an APSK modulator as in the DVB-SX standard
*/

Modulator* DVBSX_Modulator(
    const int index //!< Index identifying the modulator
    );


/*! \ingroup Interface
\brief Generate an non uniform QAM modulator as from the DVB-T standard
*/

Modulator* Nonuniform_QAM_Modulator(const int m, const int a);

/*! \ingroup Interface
\brief Generate a Modulator as specified in the xDSL standard

*/

Modulator* xDSL_Modulator(const int m);


/* \ingroup Interface
\brief Generate a Gaussian modulator.

The points of the constellation are the centroids of a Gaussian distribution.

*/

Modulator* Gaussian_Modulator(
    const int m, //!< Number of modulation bits
    const int* mapping = 0 //!< Permutation for the mapping
    );

/* \ingroup Interface
\brief Read the constellation points from a file.

The points of the constellation are read from the file. The mapping is assumed to be  natural.

*/

Modulator* Read_from_File(char* namefile);


static const double c128135180_remap[256] ={
	0.338300,0.071000,
	1.249500,0.426700,
	1.317200,0.091300,
	1.292500,0.269700,
	0.586200,0.089100,
	0.708000,0.186900,
	0.947500,0.078700,
	0.894200,0.243800,
	0.289600,0.188800,
	1.185200,0.581700,
	0.995900,0.866800,
	1.105100,0.722500,
	0.477700,0.351200,
	0.633000,0.368100,
	0.726000,0.613800,
	0.804400,0.460400,
	0.071200,0.338300,
	0.426000,1.249700,
	0.092100,1.317100,
	0.269700,1.292500,
	0.089100,0.586200,
	0.186900,0.708000,
	0.078700,0.947500,
	0.243200,0.894400,
	0.188800,0.289600,
	0.581000,1.185600,
	0.866200,0.996500,
	0.723200,1.104600,
	0.351200,0.477700,
	0.368900,0.632500,
	0.614300,0.725600,
	0.460400,0.804400,
	0.338300,-0.071000,
	1.249500,-0.426700,
	1.317200,-0.091300,
	1.292500,-0.269700,
	0.586200,-0.089100,
	0.708000,-0.186900,
	0.947500,-0.078700,
	0.894200,-0.243800,
	0.289600,-0.188800,
	1.185200,-0.581700,
	0.995900,-0.866800,
	1.105100,-0.722500,
	0.477700,-0.351200,
	0.633000,-0.368100,
	0.726000,-0.613800,
	0.804400,-0.460400,
	0.071200,-0.338300,
	0.426000,-1.249700,
	0.092100,-1.317100,
	0.269700,-1.292500,
	0.089100,-0.586200,
	0.186900,-0.708000,
	0.078700,-0.947500,
	0.243200,-0.894400,
	0.188800,-0.289600,
	0.581000,-1.185600,
	0.866200,-0.996500,
	0.723200,-1.104600,
	0.351200,-0.477700,
	0.368900,-0.632500,
	0.614300,-0.725600,
	0.460400,-0.804400,
	-1.249500,-0.426700,
	-0.338300,-0.071000,
	-0.586200,-0.089100,
	-0.708000,-0.186900,
	-1.317200,-0.091300,
	-1.292500,-0.269700,
	-0.894200,-0.243800,
	-0.947500,-0.078700,
	-1.185200,-0.581700,
	-0.289600,-0.188800,
	-0.477700,-0.351200,
	-0.633000,-0.368100,
	-0.995900,-0.866800,
	-1.105100,-0.722500,
	-0.804400,-0.460400,
	-0.726000,-0.613800,
	-0.426000,-1.249700,
	-0.071200,-0.338300,
	-0.089100,-0.586200,
	-0.186900,-0.708000,
	-0.092100,-1.317100,
	-0.269700,-1.292500,
	-0.243200,-0.894400,
	-0.078700,-0.947500,
	-0.581000,-1.185600,
	-0.188800,-0.289600,
	-0.351200,-0.477700,
	-0.368900,-0.632500,
	-0.866200,-0.996500,
	-0.723200,-1.104600,
	-0.460400,-0.804400,
	-0.614300,-0.725600,
	-1.249500,0.426700,
	-0.338300,0.071000,
	-0.586200,0.089100,
	-0.708000,0.186900,
	-1.317200,0.091300,
	-1.292500,0.269700,
	-0.894200,0.243800,
	-0.947500,0.078700,
	-1.185200,0.581700,
	-0.289600,0.188800,
	-0.477700,0.351200,
	-0.633000,0.368100,
	-0.995900,0.866800,
	-1.105100,0.722500,
	-0.804400,0.460400,
	-0.726000,0.613800,
	-0.426000,1.249700,
	-0.071200,0.338300,
	-0.089100,0.586200,
	-0.186900,0.708000,
	-0.092100,1.317100,
	-0.269700,1.292500,
	-0.243200,0.894400,
	-0.078700,0.947500,
	-0.581000,1.185600,
	-0.188800,0.289600,
	-0.351200,0.477700,
	-0.368900,0.632500,
	-0.866200,0.996500,
	-0.723200,1.104600,
	-0.460400,0.804400,
	-0.614300,0.725600};

