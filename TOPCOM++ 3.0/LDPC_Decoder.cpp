// LDPC_Decoder1.cpp: implementation of the LDPC_Decoder class.
//
//////////////////////////////////////////////////////////////////////
#include "LDPC_Decoder.h"
//#include "LDPC_Slope_Design.h"
#include "Tool_functions.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const int minrel = 0;
LDPC_Decoder::LDPC_Decoder()
{
	EXT = 0;
	maxtab = 0;
	ic = 0;
	nconnvar = 0;
	nconncheck = 0;
	perm = 0;
	buf = 0;
	gamma = 0;
	OLD = 0;//Storage of old messages for damping and/or RBP
	damp = 0;
	sat = 0;// Saturation of messages
	Lambda = 0;
	stop = true;
	checkmax = 10000000;
}
LDPC_Decoder::~LDPC_Decoder()
{
	delete[] ic;
	delete[] nconnvar;
	delete[] nconncheck;
	delete[] perm;
	delete[] buf;
	delete[] EXT;
	delete[] OLD;
	delete[] Lambda;
}
void LDPC_Decoder::Set_Precision(const double fact)
{
	delete[] ic;
	int i;
	maxtab = -1;
	int temp;
	do
	{
		maxtab++;
		temp = (int)((double)fact*log(1. + exp(-(double)maxtab / fact)) + 0.5);
	} while (temp != 0);

	if (maxtab > 0)
	{
		ic = new int[maxtab];
		for (i = 0; i < maxtab; i++)
		{
			ic[i] = (int)(fact*log(1. + exp(-i / fact)) + 0.5);
		}
	}
	else
	{
		ic = 0;
	}
}
void LDPC_Decoder::SetRBP(const double* gamma)
{
	delete[] OLD;
	this->gamma = gamma;
	if (gamma != 0)
		OLD = new int[N1];
}
void LDPC_Decoder::SetDamping(const double damp)
{
	delete[] OLD;
	this->damp = damp;
	if (damp != 0)OLD = new int[nones];
}
void LDPC_Decoder::SetParameters(const LDPC_Encoder* a, const int niterin, const double fact)
{

	int i, j, l, m;
	encoder = a;
	K = a->K;
	N1 = a->N1;	// Before puncturing
//	N =a->N;
	niter = niterin;
	int K1 = a->K1;

	delete[]nconnvar;  nconnvar = new int[N1];
	delete[]nconncheck; nconncheck = new int[N1 - K1];
	for (i = 0; i < N1; i++)nconnvar[i] = 0;
	for (i = 0; i < N1 - K1; i++)nconncheck[i] = 0;


	/* Compute variable and check node degrees */
	nones = 0;
	for (i = 0; i < (N1 - K1)*a->ncheck; i++)
	{
		if (a->Hconn[i] >= 0)
		{
			nconnvar[a->Hconn[i]]++;
			nconncheck[i / a->ncheck]++;
			nones++;
		}
		else
		{
			i = (i / a->ncheck)*a->ncheck + a->ncheck - 1;
		}
	}

	/* Compute connections of LDPC decoder */
	delete[] Lambda;
	Lambda = new int[2*N1];// For layered decoding double buffer for updating at the  next call
	delete[] perm;
	perm = new int[nones];
	int** index = new int*[N1];

	index[0] = perm;
	for (j = 1; j < N1; j++)
	{
		index[j] = index[j - 1] + nconnvar[j - 1];
		nconnvar[j - 1] = 0;
	}
	nconnvar[N1 - 1] = 0;

	l = 0;
	for (i = 0, m = 0; i < (N1 - K1)*a->ncheck; i++)// lines
	{
		if (a->Hconn[i] >= 0)
		{
			index[a->Hconn[i]][nconnvar[a->Hconn[i]]++] = m++;
		}
		else
		{
			i = (i / a->ncheck)*a->ncheck + a->ncheck - 1;
		}
	}
	delete[]   index;
	//for(i=0;i<N1;i+=2025)
	//{
	//	if(i%8100==0)printf("\n");printf("%d ",nconnvar[i]);
	//}
	maxdeg = 0;
	for (i = 0; i<N1 - K1; i++)if (nconncheck[i]>maxdeg)maxdeg = nconncheck[i];
	if (maxdeg > 1)buf = new int[maxdeg - 1];

	if (EXT != 0)delete[]  EXT;
	EXT = new int[nones];

	Set_Precision(fact);
	// Reduce number of checked parities for stopping (only for non genie aided stop)
	if (a->ccore != 0)checkmax = (a->ccore)*(a->z);
	else              checkmax = 1000000;
}
double LDPC_Decoder::Run(const int ntics, const int* LLRCh, int* out, const int* data, int*out2, bool reset)
{
	int i, ii, j, it;
	int temp, temp1, tt;
	const int* pp;
	int* ext;
	int aviter = 0;
	int jj;
	int K1 = K + encoder->shorten;
	for (ii = 0; ii < ntics; ii++)
	{
		int err = 1;
		int err1 = 10000;
		const unsigned int a = 0x80000000;
		/* initialization */
		if (reset)
			for (i = 0; i < nones; i++)
			{
				EXT[i] = 0;
			}
		for (it = 1; it <= niter; it++)
		{

			/* Variable nodes */
			pp = perm;
			err = 0;
			jj = 0;
			for (j = 0; j < encoder->shorten; j++)// Shortened bits
			{
				for (i = 0; i < nconnvar[j]; i++)EXT[pp[i]] = -1000000;
				pp += nconnvar[j];
			}
			for (; j < N1; j++)
			{
				temp = EXT[pp[0]];
				for (i = 1; i < nconnvar[j]; i++)temp += EXT[pp[i]];

				if (encoder->period == 0)
				{
					if (out2 != 0)out2[j] = temp; // Extrinsic
					temp += LLRCh[j];
					//					if(out2!=0)out2[j]=temp;

				}
				else
				{
					if (encoder->pattern[j%encoder->period])
					{
						if (out2 != 0)out2[jj] = temp;
						temp += LLRCh[jj++];
						//					if(out2!=0)out2[j]=temp;
					}
				}

				/* RBP */
				if (gamma != 0)
				{
					temp1 = temp;
					if (it > 1)temp += (int)(gamma[it - 1] * OLD[j]);
					OLD[j] = temp1;
				}

				/* Output */
				if (j < K1)
				{
					out[j - encoder->shorten] = temp;
					if (data != 0)
					{
						if (out[j - encoder->shorten] >= 0 && data[j] == 0)err++;
						else if (out[j - encoder->shorten] <= 0 && data[j] == 1)err++;
					}
					else err = 1;
				}

				for (i = 0; i < nconnvar[j]; i++)
				{
					EXT[pp[i]] = temp - EXT[pp[i]];
				}

				pp += nconnvar[j];
			}
			// Saturation
			if (sat > 0)
			{
				for (j = 0; j < N1; j++)
				{
					if (EXT[j] > +sat)EXT[j] = +sat;
					else if (EXT[j] < -sat)EXT[j] = -sat;
				}
			}
			/* Check nodes    */
			ext = EXT;
			for (j = 0; j < N1 - K1; j++)
			{
				if (nconncheck[j] > 2)
				{
					buf[nconncheck[j] - 2] = ext[nconncheck[j] - 1];
					for (i = nconncheck[j] - 3; i >= 0; i--)
					{
						buf[i] = g(buf[i + 1], ext[i + 1]);
					}

					tt = ext[0];
					ext[0] = buf[0];
					if (data == 0)	if ((ext[0] ^ tt) <= 0)err++; // Check sign difference


					for (i = 1; i < nconncheck[j] - 1; i++)
					{
						temp = ext[i];
						ext[i] = g(tt, buf[i]);
						if (data == 0)if ((ext[i] ^ temp) < 0)err++; // Check sign
						tt = g(tt, temp);
					};
					if (data == 0)if ((ext[nconncheck[j] - 1] ^ tt) < 0)err++;
					ext[nconncheck[j] - 1] = tt;
				}
				else if (nconncheck[j] == 2)
				{
					if (data == 0)if ((((ext[0] ^ ext[1]) & a) != 0))err++; // Check sign difference
					tt = ext[0];
					ext[0] = ext[1];
					ext[1] = tt;
				}
				else if (nconncheck[j] == 1)
				{
					ext[0] = -100000;
				}
				ext += nconncheck[j];
			}
			//			if(err<12)break;
			if (err == 0)break;
			//if(err>err1)
			//  printf("%d %d %d\n",it,err1,err);
			//err1=err;
		}
		LLRCh += encoder->N;
		out += K;
		if (out2 != 0)out2 += encoder->N;
		if (data != 0)data += K;
		reset = true;
		aviter += it;
	}
	return (double)aviter / ntics;
}
double LDPC_Decoder::RunSubopt(const int ntics, const int* LLRCh, int* out, const int* data, const int beta, int*out2, bool reset)
{
	int i, ii, j, it, sign, min, minabs, min2abs, ab;
	int temp, temp1;
	const int* pp;
	int* ext;
	int aviter;
	int jj;
	int err;
	aviter = 0;
	int K1 = K + encoder->shorten;
	for (ii = 0; ii < ntics; ii++)
	{

		err = 1;
		/* Initialization */
		if (reset)
			for (i = 0; i < nones; i++)EXT[i] = 0;


		for (it = 1; it <= niter; it++)
		{
			pp = perm;
			err = 0;
			jj = 0;

			for (j = 0; j < encoder->shorten; j++)// Shortened bits
			{
				for (i = 0; i < nconnvar[j]; i++)EXT[pp[i]] = -1000000;
				pp += nconnvar[j];
			}
			/* Variable nodes */
			for (; j < N1; j++)
			{
				temp = EXT[pp[0]];
				for (i = 1; i < nconnvar[j]; i++)temp += EXT[pp[i]];
				if (encoder->period == 0)
				{
					temp += LLRCh[jj];
					if (out2 != 0)out2[jj] = temp;
					jj++;
				}
				else
				{
					if (encoder->pattern[j%encoder->period])
					{
						if (out2 != 0)out2[jj] = temp;// Extrinsic;
						temp += LLRCh[jj++];
					}
				}
				if (j < K1)
				{
					out[j - encoder->shorten] = temp;
					if (data != 0)
					{
						if (out[j - encoder->shorten] >= 0 && data[j] == 0)err++;
						else if (out[j - encoder->shorten] <= 0 && data[j] == 1)err++;
					}
					else err = 1;
				}
				/* RBP */
				if (gamma != 0)
				{
					temp1 = temp;
					temp += (int)(gamma[it] * OLD[j]);
					OLD[j] = temp1;
				}

				// Update extrinsic...
				for (i = 0; i < nconnvar[j]; i++)EXT[pp[i]] = temp - EXT[pp[i]];
				pp += nconnvar[j];
			}


			/* Check nodes    */
			ext = EXT;
			for (j = 0; j < N1 - K1; j++)
			{
				if (nconncheck[j] == 1)ext[0] = -100000;
				else
				{
					sign = 0;
					minabs = min2abs = 1000000000;
					for (i = 0; i < nconncheck[j]; i++)
					{
						sign ^= -ext[i];
						if (ext[i] > 0)ab = ext[i];
						else	    ab = -ext[i];
						if (ab < minabs)
						{
							min2abs = minabs;
							minabs = ab;
							min = i;
						}
						else if (ab < min2abs)
						{
							min2abs = ab;
						}
					}
					// Offset correction
					if (minabs > beta) minabs -= beta;  else minabs = 0;
					if (min2abs > beta)min2abs -= beta; else min2abs = 0;
					//					printf("%d %d,%d\n",min,minabs,min2abs);
					for (i = 0; i < nconncheck[j]; i++)
					{
						if (i == min)
						{
							if ((sign ^ (-ext[i])) >= 0)ext[i] = -min2abs;
							else			       ext[i] = min2abs;
						}
						else
						{
							if ((sign ^ (-ext[i])) >= 0) ext[i] = -minabs;
							else					ext[i] = minabs;
						}
					}
				}
				ext += nconncheck[j];
			}

			if (err == 0)break;
		}
		LLRCh += encoder->N;
		out += K;
		if (out2 != 0)out2 += encoder->N;
		if (data != 0)data += K;
		reset = true;
		aviter += it;
	}
	return (double)aviter / ntics;
}
double LDPC_Decoder::RunDamped(const int ntics, const int* LLRCh, int* out, const int* data, const int beta, int*out2, bool reset)
{
	int i, ii, j, it, sign, min, minabs, min2abs, ab;
	int temp;
	const int* pp;
	int* ext;
	int aviter;
	int jj;
	int err;
	aviter = 0;
	for (ii = 0; ii < ntics; ii++)
	{

		err = 1;
		const unsigned int a = 0x80000000;

		/* Initialization */
		double tempR;

		for (it = 1; it <= niter; it++)
		{
			pp = perm;
			err = 0;
			jj = 0;

			/* Variable nodes */
			for (j = 0; j < N1; j++)
			{
				temp = 0;
				if (!reset)
				{
					tempR = (1. - damp)*(double)OLD[pp[0]] + damp*(double)EXT[pp[0]];
					for (i = 1; i < nconnvar[j]; i++)
					{
						tempR += (1. - damp)*(double)OLD[pp[i]] + damp*(double)EXT[pp[i]];
					}
					if (tempR > 0)	temp = (int)(tempR + 0.5);
					else		temp = (int)(tempR - 0.5);
				}

				/* Add channel values  */
				if (encoder->period == 0)
				{
					if (out2 != 0)out2[j] = temp;
					temp += LLRCh[j];
				}
				else
				{
					if (encoder->pattern[j%encoder->period])
					{
						if (out2 != 0)out2[jj] = temp;
						temp += LLRCh[jj++];
					}
				}


				if (j < K)
				{
					out[j] = temp;
					if (data != 0)
					{
						if (out[j] >= 0 && data[j] == 0)
						{
							err++;
						}
						else if (out[j] <= 0 && data[j] == 1)
						{
							err++;
						}
					}
					else err = 1;
				}

				// Update extrinsic...
				if (reset)
				{
					for (i = 0; i < nconnvar[j]; i++)
					{
						OLD[pp[i]] = 0;
						EXT[pp[i]] = temp;
					}
				}
				else
				{
					for (i = 0; i < nconnvar[j]; i++)
					{
						OLD[pp[i]] = EXT[pp[i]];
						EXT[pp[i]] = temp - EXT[pp[i]];
					}
				}
				pp += nconnvar[j];
			}
			reset = false;



			/* Check nodes    */
			ext = EXT;
			for (j = 0; j < N1 - K; j++)
			{
				if (nconncheck[j] == 1)ext[0] = -100000;
				else
				{
					sign = 0;
					minabs = min2abs = 1000000000;
					for (i = 0; i < nconncheck[j]; i++)
					{
						sign ^= -ext[i];
						if (ext[i] > 0)ab = ext[i];
						else	    ab = -ext[i];
						if (ab < minabs)
						{
							min2abs = minabs;
							minabs = ab;
							min = i;
						}
						else if (ab < min2abs)
						{
							min2abs = ab;
						}
					}
					// Offset correction
					if (minabs > beta) minabs -= beta;  else minabs = 0;
					if (min2abs > beta)min2abs -= beta; else min2abs = 0;
					//					printf("%d %d,%d\n",min,minabs,min2abs);
					for (i = 0; i < nconncheck[j]; i++)
					{
						if (i == min)
						{
							if ((sign ^ (-ext[i])) >= 0)ext[i] = -min2abs;
							else			       ext[i] = min2abs;
						}
						else
						{
							if ((sign ^ (-ext[i])) >= 0) ext[i] = -minabs;
							else					ext[i] = minabs;
						}
					}
				}
				ext += nconncheck[j];
			}
			if (err == 0)break;
		}
		LLRCh += encoder->N;
		out += K;
		if (out2 != 0)out2 += encoder->N;
		if (data != 0)data += K;
		reset = true;
		aviter += it;
	}
	return (double)aviter / ntics;
}
double LDPC_Decoder::RunLayered(const int ntics, const int* LLRCh, int* out, const int* data, bool opt, const int beta, int*out2, bool reset)
{
	int tt, temp,tt1;
	int i, ii, ll, j, it, sign, sign1, min, minabs, min2abs, ab;
	const int* pp;
	int* ext;
	int aviter;
	int jj;
	int err, t;
	aviter = 0;
	int K1 = K + encoder->shorten;
	for (ii = 0; ii < ntics; ii++)
	{
		if (reset)// Initialize messages and load Channel values with possible shortening and puncturing
		{
			for (i = 0; i < nones; i++)EXT[i] = 0;
			for (i = jj=ll=0; i < N1; i++)
			{
				if (encoder->lastsh == 0 && i < encoder->shorten)Lambda[i] = -100000;
				else if (encoder->lastsh == 1 && (i >= K && i < K + encoder->shorten))Lambda[i] = -100000;
				else // Puncturing pattern is applied *after* shortening
				{
					if (encoder->period == 0)Lambda[i] = LLRCh[jj++];
					else
					{
						if (encoder->pattern[ll%encoder->period])Lambda[i] = LLRCh[jj++];
						else								     Lambda[i] = 0;
						ll++;
					}

				}
				Lambda[i + N1] = Lambda[i];
			}
		}
		else
		{
			for (i = jj=ll=0; i < N1; i++)
			{
				if (encoder->lastsh == 0 && i < encoder->shorten)Lambda[i] = -100000;
				else if (encoder->lastsh == 1 && (i >= K && i < K + encoder->shorten))Lambda[i] = -100000;
				else // Puncturing pattern is applied *after* shortening
				{
					if (encoder->period == 0)
					{
						Lambda[i]     += (LLRCh[jj] - Lambda[i + N1]);
						Lambda[i + N1] = (LLRCh[jj++]);
					}

					else
					{
						if (encoder->pattern[ll%encoder->period])
						{
							Lambda[i] += (LLRCh[jj] - Lambda[i + N1]);
							Lambda[i + N1] = (LLRCh[jj++]);
						}
						else								     Lambda[i] = 0;
						ll++;
					}
				}
			}
		}

		for (it = 1; it <= niter; it++)
		{
			err = 0;
			ext = EXT;
			if (!opt) // Min-sum-offset
			{
				for (j = 0; j < N1 - K1; j++)// Check Node processing
				{
					pp = encoder->Hconn + j*encoder->ncheck;
					minabs = min2abs = 1000000000;
					sign = 0;
					for (i = 0; i < nconncheck[j]; i++)
					{
						t = (Lambda[pp[i]] -= ext[i]); // Temporary value
						if (t > 0) { ab = t; sign ^= 1; }
						else { ab = -t; }
						if (ab < minabs)
						{
							min2abs = minabs;
							minabs = ab;
							min = i;
						}
						else if (ab < min2abs)min2abs = ab;
					}
					if (data == 0 && j < checkmax)
					{
						if(sign == 1 )err++;
						else if(minabs<=minrel)err++;
					}
					// Offset correction
					if (minabs > beta) minabs -= beta;  else minabs = 0;
					if (min2abs > beta)min2abs -= beta; else min2abs = 0;
					if (sat > 0)
					{
						if (minabs > sat)	minabs = sat;
						if (min2abs > sat)	min2abs = sat;
					}

					for (i = 0; i < nconncheck[j]; i++)
					{
						if (Lambda[pp[i]]>0) sign1 = sign ^ 1;
						else				sign1 = sign;
						if (i == min)
						{
							if (sign1)ext[i] = min2abs;
							else	 ext[i] = -min2abs;
						}
						else
						{
							if (sign1)ext[i] = minabs;
							else	 ext[i] = -minabs;
						}
						Lambda[pp[i]] += ext[i];

					}
					ext += nconncheck[j];
				}
			}
			else     // Optimal processing
			{
				for (j = 0; j < N1 - K1; j++)// Check Node processing
				{
					pp = encoder->Hconn + j*encoder->ncheck;
					if (nconncheck[j] > 2)
					{
						buf[nconncheck[j] - 2] = (Lambda[pp[nconncheck[j] - 1]] -= ext[nconncheck[j] - 1]); // //= ext[nconncheck[j]-1];
						for (i = nconncheck[j] - 3; i >= 0; i--)
						{
							buf[i] = g(buf[i + 1], Lambda[pp[i + 1]] - ext[i + 1]);
						}
						if (sat == 0)
						{
							tt = (Lambda[pp[0]] -= ext[0]);
							ext[0] = buf[0];
							Lambda[pp[0]] += ext[0];
							for (i = 1; i < nconncheck[j] - 1; i++)
							{
								temp = (Lambda[pp[i]] -= ext[i]);//ext[i];
								ext[i] = g(tt, buf[i]);
								Lambda[pp[i]] += ext[i];
								tt = g(tt, temp);
							}
							ext[nconncheck[j] - 1] = tt;
							if (data == 0 && j<checkmax)
							{
								tt1 = g(tt,Lambda[pp[nconncheck[j] - 1]]);
//								if (tt1 >= 0)
								if (tt1 >= -minrel)err++;

							}
							Lambda[pp[nconncheck[j] - 1]] += tt;
						}
						else
						{
							tt = (Lambda[pp[0]] -= ext[0]);
							ext[0] = buf[0];
							if (ext[0] > sat)ext[0] = +sat;
							else if (ext[0] < -sat)ext[0] = -sat;
							Lambda[pp[0]] += ext[0];
							for (i = 1; i < nconncheck[j] - 1; i++)
							{
								temp = (Lambda[pp[i]] -= ext[i]);//ext[i];
								ext[i] = g(tt, buf[i]);
								if (ext[i] > sat)ext[i] = +sat;
								else if (ext[i]<-sat)ext[i] = -sat;
								Lambda[pp[i]] += ext[i];
								tt = g(tt, temp);
							};
							ext[nconncheck[j] - 1] = tt;
							if (ext[nconncheck[j] - 1]> sat)ext[nconncheck[j] - 1] = +sat;
							else if (ext[nconncheck[j] - 1] < -sat)ext[nconncheck[j] - 1] = -sat;
							if (data == 0 && j < checkmax)
							{
								tt = g(tt, Lambda[pp[nconncheck[j] - 1]]);
//								if (tt >= 0)err++;
								if (tt1 >= -minrel)err++;
							}
							Lambda[pp[nconncheck[j] - 1]] += ext[nconncheck[j] - 1];
						}
					}
					else if (nconncheck[j] == 2)
					{
						tt = Lambda[pp[0]] + Lambda[pp[1]] - ext[0] - ext[1];

						t    = Lambda[pp[1]] - ext[1];
						tt1  = Lambda[pp[0]] - ext[0];
						ext[1] = tt1;
						ext[0] = t;

						Lambda[pp[0]] = Lambda[pp[1]] = tt;
						if (data == 0 && j<checkmax)
						{
							tt1 = g(tt1, t);
							if (tt1 >= -minrel)err++;
						}
					}
					else if (nconncheck[j] == 1)
					{
						ext[0] = -100000;
						Lambda[pp[0]] = -100000;
					}
					ext += nconncheck[j];
				}
			}
			// Write Output and check is valid
			for (j = 0; j < K; j++)
			{
				if (encoder->lastsh == 0)out[j] = Lambda[j + encoder->shorten];
				else                     out[j] = Lambda[j];
				if (data != 0)
				{
					if ((out[j] >= 0 && data[j] == 0) || (out[j] <= 0 && data[j] == 1))
					{
//						printf("%d->%d\n", it, j);
						err++;
					}

				}
			}
//			printf("it=%d\terr=%d\n", it, err);
			if(err == 0 && stop)break;
		}

		if (out2 != 0)// Store LLR of coded bit in optional buffer
		{
			for (i = jj = ll = 0; i < N1; i++)
			{
				if (!((encoder->lastsh == 0 && i < encoder->shorten)
						||(encoder->lastsh == 1 && (i >= K && i < K + encoder->shorten))
							))
				{
					if (encoder->period == 0)out2[jj++]=Lambda[i]-Lambda[i+N1];  // Extrinsic
					else
					{
						if (encoder->pattern[ll%encoder->period])out2[jj++] = Lambda[i]-Lambda[i+N1];
						ll++;
					}
				}
			}
			out2 += encoder->N;
		}
		LLRCh += encoder->N;
		out += K;
		if (data != 0)data += K;
		reset = true;
		aviter += it;
	}
	return (double)aviter / ntics;
}
double LDPC_Decoder::RunLayeredQ(const int ntics, const int* LLRCh, int* out, const int* data, bool opt, const int beta)
{
	int tt;
	int i, ii, j, it, sign, sign1, min, minabs, min2abs, ab;
	const int* pp;
	int* ext;
	int aviter;
	int jj;
	const int NLLR = 3;
	int err, t;
	aviter = 0;
	int K1 = K + encoder->shorten;

	for (ii = 0; ii < ntics; ii++)
	{
		for (i = 0; i < nones; i++)EXT[i] = 0; // accumulators
		for (i = jj = 0; i < N1 / 2; i++)
		{
			if (encoder->period == 0)
			{
				for (j = 0; j < NLLR; j++)Lambda[i*NLLR + j] = LLRCh[jj++];
			}

			else
			{
				if (encoder->pattern[i%encoder->period])
				{
					for (j = 0; j < NLLR; j++)Lambda[i*NLLR + j] = LLRCh[jj++];
				}
				else
				{
					for (j = 0; j< NLLR; j++)Lambda[i*NLLR + j] = 0;
				}
			}
		}
		int* ppp;
		int lll[60];
		int add;
		for (it = 1; it <= niter; it++)
		{
			err = 0;
			ext = EXT;
			if (!opt) // Min-sum-offset
			{
				for (j = 0; j < N1 - K1; j++)// Check Node processing
				{
					pp = encoder->Hconn + j*encoder->ncheck;
					minabs = min2abs = 1000000000;
					sign = 0;
					for (i = 0; i < nconncheck[j]; i++)
					{
						add = encoder->permvar[pp[i]];// psym[pp[i]];// map to symbol bit
						ppp = Lambda + 3 * (add / 2);
						if (add & 1)//odd
						{
							t = lll[i]= maxx(ppp[1], ppp[2]) - maxx(ppp[0], 0) - ext[i];
							ppp[1] -= ext[i];
							ppp[2] -= ext[i];
						}
						else//even
						{
							t  = lll[i]= maxx(ppp[0], ppp[2]) - maxx(ppp[1], 0) - ext[i];
							ppp[0] -= ext[i];
							ppp[2] -= ext[i];
						}
						//						t = (Lambda[pp[i]] -= ext[i]); // Temporary value
						if (t > 0) { ab = t; sign ^= 1; }
						else { ab = -t; }
						if (ab < minabs)
						{
							min2abs = minabs;
							minabs = ab;
							min = i;
						}
						else if (ab < min2abs)min2abs = ab;
					}
					if (data == 0 && sign == 1 && j<checkmax)err++;
					// Offset correction
					if (minabs > beta) minabs -= beta;  else minabs = 0;
					if (min2abs > beta)min2abs -= beta; else min2abs = 0;
					if (sat > 0)
					{
						if (minabs > sat)	minabs = sat;
						if (min2abs > sat)	min2abs = sat;
					}

					for (i = 0; i < nconncheck[j]; i++)
					{
						if (lll[i]>0) sign1 = sign ^ 1;
						else		   sign1 = sign;
						if (i == min)
						{
							if (sign1)ext[i] = min2abs;
							else	 ext[i] = -min2abs;
						}
						else
						{
							if (sign1)ext[i] = minabs;
							else	 ext[i] = -minabs;
						}
						add = encoder->permvar[pp[i]];
						ppp = Lambda + 3 * (add / 2);
						if (add & 1)//odd
						{
							ppp[1] += ext[i];
							ppp[2] += ext[i];
						}
						else//even
						{
							ppp[0] += ext[i];
							ppp[2] += ext[i];
						}
					}
					ext += nconncheck[j];
				}
			}
			else     // Optimal processing
			{
				int deg;
				for (j = 0; j < N1 - K1; j++)// Check Node processing
				{
					pp  = encoder->Hconn + j*encoder->ncheck;
					deg = nconncheck[j];
					if (deg <= 2)printf("Warning: %d\n", deg);
					for (i = deg - 1; i >= 0; i--)
					{
						add = encoder->permvar[pp[i]];
						ppp = Lambda + 3 * (add / 2);
						if (add & 1)//odd
						{
							lll[i] = maxx(ppp[1], ppp[2]) - maxx(ppp[0], 0) - ext[i];
							ppp[1] -= ext[i];
							ppp[2] -= ext[i];
						}
						else//even
						{
							lll[i] = maxx(ppp[0], ppp[2]) - maxx(ppp[1], 0) - ext[i];
							ppp[0] -= ext[i];
							ppp[2] -= ext[i];
						}
						if(i==deg-1)buf[i-1] = lll[i];
						else if(i>0)buf[i-1] = g(buf[i], lll[i]);
					}
					for (i = 0; i < deg ; i++)
					{
						if(i==0)ext[i] = buf[0];
						else if (i == deg - 1)ext[i] = tt;
						else   ext[i] = g(tt, buf[i]);
						add = encoder->permvar[pp[i]];
						ppp = Lambda + 3 * (add / 2);
						if (add & 1)//odd
						{
							ppp[1] += ext[i];
							ppp[2] += ext[i];
						}
						else//even
						{
							ppp[0] += ext[i];
							ppp[2] += ext[i];
						}							
						if(i==0)tt = lll[i];
						else    tt = g(tt, lll[i]);
					}
					if (data == 0 && j<checkmax)
					{
						if (tt >= 0)err++;
					}
					ext += nconncheck[j];
				}
			}
			// Write Output and check is valid
			for (j = 0; j < K; j++)
			{
				if (encoder->lastsh == 0)out[j] = Lambda[j + encoder->shorten];
				else                     out[j] = Lambda[j];
				if (data != 0)
				{
					if ((out[j] >= 0 && data[j] == 0) || (out[j] <= 0 && data[j] == 1))
					{
						//						printf("%d->%d\n", it, j);
						err++;
					}

				}
			}
			//			printf("it=%d\terr=%d\n", it, err);
			if (err == 0 && stop)break;
		}
		LLRCh += encoder->N;
		out += K;
		if (data != 0)data += K;
		aviter += it;
	}
	return (double)aviter / ntics;
}
//#ifdef WIN32
//__forceinline int LDPC_Decoder::g(const int a, const int b)
//#else
void LDPC_Decoder::Display(FILE* file) const
{
	int i;
	int *stat = (int*)calloc(1000, sizeof(int));
	double tot = 0;
	int K1 = K + encoder->shorten;
	for (i = 0; i < N1; i++)
	{
		stat[nconnvar[i]]++;
		tot += nconnvar[i];
	}
	for (i = 0; i < 1000; i++)
	{
		if (stat[i] == 0)continue;
		fprintf_s(file, "var.\t%2d\t%d \t %f\n", i, stat[i], (double)stat[i] / N1);
		stat[i] = 0;
	}
	fprintf_s(file, "\n");
	for (i = 0; i < N1 - K1; i++)
	{
		stat[nconncheck[i]]++;
		tot += nconncheck[i];
	}
	for (i = 0; i < 1000; i++)
	{
		if (stat[i] == 0)continue;
		fprintf_s(file, "chk.\t%2d\t%d \t%f\n", i, stat[i], (double)stat[i] / (N1 - K));
		stat[i] = 0;
	}
	fprintf_s(file, "avev=\t%f\t avec=\t%f\n", (double)nones / N1, (double)nones / (N1 - K));
	delete[] stat;
}

void LDPC_Decoder::LongDisplay(FILE* file) const
{
	int i,ll;
	int Z = encoder->GetZ();
	int K1 = K + encoder->K1;
	fprintf(file, "VN_________________________\n");
	for (i = ll=0; i < N1; i++)
	{
//		fprintf(file, "%d\t%d\t", i, nconnvar[i]);
		if (encoder->lastsh == 0 && i < encoder->shorten)fprintf(file,"S");
		else if (encoder->lastsh == 1 && (i >= K && i < K + encoder->shorten))fprintf(file, "S");
		else // Puncturing pattern is applied *after* shortening
		{
			if (encoder->period !=0)
			{
				if (encoder->pattern[ll%encoder->period]!=1)fprintf(file,"P");					     
				else fprintf(file, "X");		
				ll++;
			}
			else fprintf(file, "X");		
		}
		if(i%Z==Z-1)fprintf(file, "\n");
	}
	fprintf(file, "CN__________________________\n");
	int j;
	for (i = 0; i < N1 - K1; i+=Z)
	{
		fprintf(file, "%d\t%d\t", i,nconncheck[i]);
		for(j=0;j<nconncheck[i];j++)
		fprintf(file, "%d\t%d\t", encoder->Hconn[i*encoder->ncheck+j]/Z, encoder->Hconn[i*encoder->ncheck + j] % Z);
		fprintf(file, "\n");
	}
}


inline int LDPC_Decoder::maxx(const int a, const int b)
{
	int t, c;
	t = (a - b);
	if (t > 0) { c = a; }
	else { c = b; t = -t; }
	if (t < maxtab)c += ic[t];
	return c;
}
inline int LDPC_Decoder::g(const int a, const int b)
//#endif
{
	int t, c;

	t = (a - b);
	if (t > 0) { c = a; }
	else { c = b; t = -t; }
	if (t < maxtab)c += ic[t];

	t = (a + b);
	if (t > 0) { c -= t; }
	else { t = -t; }
	if (t < maxtab)c -= ic[t];

	return c;
}
