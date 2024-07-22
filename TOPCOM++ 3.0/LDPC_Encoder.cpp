#include "LDPC_Encoder.h"
//#include "LDPC_Slope_design.h"
// LDPC_Encoder.cpp: implementation of the LDPC_Encoder class.
//
//////////////////////////////////////////////////////////////////////

#include "LDPC_Encoder.h"
#include "LDPC_Slope_Design.h"
#include "Tool_functions.h"
#include "Random_functions.h"
#include "Biggirth.h"
#include "Interleaver.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
LDPC_Encoder::LDPC_Encoder()
{
	K=K1=N=N1=0;
	counter = 0;
	G=0;
	M=2;
	shorten=0;

	period=0;
	pattern=0;
	enctype=0;
	HARQ=false;
	ccore = 0;
	z = 1;
	hhh = 0;
	aaa = 0;
	lastsh = 0;

	buff=0;
	bufHARQ=0;
	G = 0;
	Hconn = 0;	
	Gconn = 0;
	p=0;
	order=0;
	gg = 0;
	gap = 0;
	permvar = 0;
}
LDPC_Encoder::~LDPC_Encoder()
{
	delete[] Hconn;
	delete[] G;
	delete[] pattern;
	delete[] p;
	delete[] buff;
	delete[] order; 
	delete[] bufHARQ; 
	delete[] gg; gg = 0;
	delete[] permvar;

}
void LDPC_Encoder::Display(FILE* f, const int Z)const
{
	int i,j,k,l,a;
	if (Z != 0)
	{
		fprintf_s(f,"#K1=%d,N1=%d,N=%d,Z=%d, dcmax=%d\n",K1/Z, N1 / Z, N / Z, Z,ncheck);
		for (i = 0; i<(N1 - K1); i += Z)
		{
			l = 0;
			for (j = 0; j<ncheck; j++)
			{
				if (Hconn[i*ncheck + j]<0)break;
				a = Hconn[i*ncheck + j] / Z;
				for (k = l; k < a;k++)fprintf_s(f, "\t");
				fprintf_s(f, "%d\t", Hconn[i*ncheck + j] % Z);				
				l = a+1;
			}
			fprintf_s(f, "\n");
		}

	}
	else
	{
		fprintf_s(f,"#K=%d,N=%d,dcmax=%d\n",K,N,ncheck);
		for(i=0;i<(N1-K1);i++)
		{
			for(j=0;j<ncheck;j++)
			{
				fprintf_s(f,"%d\t",Hconn[i*ncheck+j]);
				if(Hconn[i*ncheck+j]<0)break;
			}
			fprintf_s(f,"\n");
		}
	}
}

void LDPC_Encoder::GetCode(int &N, int &E, int &Z, int &dC, int* &P, int* &add, int* &sh, int* &deg) const
{
	int i, j,ll;
	Z = GetZ();
	deg = new int[(N1 - K1) / Z];
	P = new int[N1/Z];
	N = N1 / Z;
	E = 0;
	dC = 0;
	for (i = 0; i < N1 - K1; i += Z)
	{
		for (j = 0; j < ncheck; j++)
		{
			if (Hconn[i*ncheck + j] < 0)break;
			E++;
		}
		deg[i / Z] = j;
		if (dC < deg[i / Z])dC = deg[i / Z];
	}

	for (i = 0; i < N1 / Z; i++)P[i] = 0;
	for (i = ll = 0; i < N1; i++)
	{
		if (lastsh == 0 && i < shorten)P[i/Z] -=1;
		else if (lastsh == 1 && (i >= K && i < K + shorten))P[i/Z] -=1;
		else // Puncturing pattern is applied *after* shortening
		{
			if (period != 0)
			{
				if (pattern[ll%period] != 1)P[i/Z] += 1;
				else P[i/Z] += 0;
				ll++;
			}
			else P[i/Z] +=0;
		}
	}

	add = new int[E];
	sh = new int[E];
	for (i = E= 0; i < N1 - K1; i += Z)
	{
		for (j = 0; j < ncheck; j++)
		{
			if (Hconn[i*ncheck + j] < 0)break;
			add[E] = Hconn[i*ncheck + j] / Z;
			sh[E] =  Hconn[i*ncheck + j] % Z;
			E++;
		}
	}
}
void LDPC_Encoder::SetParameters(const int kin, const int nin, const int ncheckin, const int *H)
{
	int i;
	K=kin;
	N=N1=nin;
	ncheck=ncheckin;
	if(Hconn!=0)delete[] Hconn;
	Hconn=new int[ncheck*(N1-K)];
	bool low_tri=true;
	for(i=0;i<(N1-K)*ncheck;i++)
	{
		Hconn[i]=H[i];
		if(Hconn[i]>K+i)low_tri=false;
	}
	if(!low_tri)// If not low triangular build generating matrix
	{
		this->SetEncodingMethod(1);
		return;  // Skip the following part
		// Build G
		M=N/32;
		if(N%32!=0)M++;
		unsigned int* HH=(unsigned int*)calloc((N-K)*M,4);
		unsigned int* h=HH;

		// Put the matrix in binary form
		for(i=0;i<(N1-K)*ncheck;i++)
		{
			if(Hconn[i]>=0)setbit(h,Hconn[i],1);
			if(i%ncheck==ncheck-1)h+=M;
		}

		// Gaussian elimination
		int* p=Gaussian_Elimination(HH,N,N-K);
		int j;

		// Printout
		for(i=0;i<N-K;i++)
		{
			for(j=0;j<N;j++)fprintf_s(stdout,"%1d",getbit(HH+i*M,p[j]));
			fprintf_s(stdout,"\n");
		}


		int* G=new int [(N1-K)*K];
		for(i=0;i<N-K;i++)
		{
			for(j=0;j<K;j++)
			{
				G[i*K+j]=getbit(HH+i*M,p[N1-K+j]);
				fprintf_s(stdout,"%1d",G[i*K+j]);	
			}
			fprintf_s(stdout,"\n");
		}
 

	}
}
void LDPC_Encoder::SetParameters(const char *namefile)
{
	delete[] Hconn; Hconn = 0;
	delete[] G; G = 0;
	delete[] pattern; pattern = 0;
	delete[] p; p = 0;
	delete[] buff; buff = 0;
	delete[] order; order = 0;
	delete[] bufHARQ; bufHARQ = 0;
	LDPC_Encoder();
	int i, j;
	int r;
	FILE* a;
	fopen_s(&a, namefile, "r");
	if (a == 0)
	{
		printf("LDPC_From_File: File %s not found.\n", namefile);
		return;
	}

	fscanf_s(a, "%d", &this->N);
	fscanf_s(a, "%d", &r);
	fscanf_s(a, "\n");
	fscanf_s(a, "%d", &this->ncheck);
	fscanf_s(a, "\n");
	this->K = this->N - r;
	this->N1 = this->N;
	this->Hconn = new int[this->N*this->ncheck];

	int*h = this->Hconn;

	int temp;
	for (i = 0; i<r; i++)
	{
		for (j = 0;; j++)
		{
			fscanf_s(a, "%d", &temp);
			if (temp<0)
			{
				h[j] = -1;
				break;
			}
			else
			{
				h[j] = temp;
			}
		}
		h += this->ncheck;
		fscanf_s(a, "\n");
	}
	fclose(a);
}
void LDPC_Encoder::Store_Matrix(const char *namefile) const
{
	int i,j;
	FILE* a;
	a=fopen(namefile, "w");
	if (a == 0)
	{
		printf("LDPC_Encoder::Store_Matrix: Cannot open %s for writing matrix.\n", namefile);
		return;
	}
	fprintf(a, "#\tN1");
	fprintf(a, "\tM1");
	fprintf(a, "\tN");
	fprintf(a, "\tK");
	fprintf(a, "\tZ");
	fprintf(a, "\tgap");
	fprintf(a, "\tmaxdc");
	fprintf(a, "\n");
	fprintf(a, "#\t%d", N1);   // Total Columns
	fprintf(a, "\t%d", N1 - K1); // Total checks
	fprintf(a, "\t%d", N);
	fprintf(a, "\t%d", K);
	fprintf(a, "\t%d", z);
	fprintf(a, "\t%d", gap);
	fprintf(a, "\t%d", ncheck);
	fprintf(a, "\n");
	int r = N1 - K1;

	int*h = this->Hconn;

	int temp;
	for (i = 0; i<r; i+=z)
	{
		for (j = 0;j<ncheck; j++)
		{
			temp =	h[j];
			fprintf(a, "%d\t", temp);
			if (temp < 0)break;
		}
		h += z*this->ncheck;
		fprintf_s(a, "\n");
	}
	if (N1>N)// Store Puncturing pattern
	{
		for (i = 0; i < N1-shorten; i++)
		{
			fprintf(a, "%1d", pattern[i%period]);
			if (i % 100 == 99)fprintf(a, "\n");
		}
		fprintf(a, "\n");
	}
	if(gg!=0)
	{ 
		int k;
		for (i = 0; i < gap; i++)
		{
			for (k = temp= 0; k < K1; k++)
			{
				temp += (gg[i*K1 + k] << (k % 32));
				if (k % 32 == 31) {fprintf(a, "%08x,", temp); temp = 0; }
			}
			fprintf(a, "\n");
		}
	}
	fclose(a);
}

// Applies after shortening
void LDPC_Encoder::SetPuncturing(const int periodin,const int *patternin)
{
	int i;
	if (patternin == 0)
	{
		period= N1 - (K1 - K);
		delete[] pattern;
		pattern=new int[period];
		for (i = 0; i < periodin; i++)pattern[i] = 0;
		for (; i < period; i++)pattern[i] = 1;
	}
	else
	{
		period=periodin;
		delete[] pattern;
		pattern=new int[period];
		for (i = 0; i < period; i++)pattern[i] = patternin[i];

	}
	N=0;
	for(i=0;i<N1-(K1-K);i++)
	{
		N+=pattern[i%period];
	}
	delete[] buff;buff=new int[N1]; // Buffer to store unpunctured bits
}
bool LDPC_Encoder::IsLowerTriangular() const
{
	int i, j;
	int*h;
	for (i = 0; i < N1 - K; i++)
	{
		h = Hconn + i*ncheck;
		for (j = 0; j < ncheck; j++)
		{
		
			if (h[j] < 0) {break; }
			if (h[j] > K + i)
			{
				return false;
			}
		}
		if (h[j-1] != K + i)return false;
	}
	return true;
}
bool LDPC_Encoder::IsCodeword(const int* inp) const
{
	int j,l,check;
	int* h;

	for(j=0;j<N1-K1;j++)
	{
		l=0;
		check=0;
		h=Hconn + ncheck*j;
		while(h[l]!=-1 && l<ncheck)
		{
//			if (h[l] >= N)goto next;// skip check if bit has been punctured  
			check+= inp[h[l]];
			l++;
		}
		check%=M;
		if(check!=0)
		{
			printf("Error:%d\n",j);
//			return false;
		}
	}
	return true;
}

void LDPC_Encoder::MakeEncodable()
{
	if (enctype == 5)return;
	int *HconnT;
	int* pcol;
	int n, n1, i, j, l, k, ncmin, nc, temp, cbest;
	if (N1%z != 0 || K1%z != 0)exit(1);
	int Nz = N1 / z;
	int Kz = K1 / z;
	int degmin;

	// Computes the left degree profile
	int* nconnvar = new int[N1/z];
	for (i = 0; i<N1/z; i++)nconnvar[i] = 0;
	int*h;
	for (i = 0; i<(N1 - K1); i+=z)
	{
		h = Hconn + i*ncheck;
		for (k = 0; k < ncheck; k++)
		{
			if (h[k] >= 0)
			{
				nconnvar[h[k]/z]++;
			}
			else break;
		}
	}
	int maxdeg = 0;
	for (i = 0; i<N1/z; i++)
	{
		if (nconnvar[i]>maxdeg)maxdeg = nconnvar[i];
	}

	// Construct the compressed matrix of columns (HconnT)
	HconnT = new int[maxdeg*Nz];
	for (i = 0; i < maxdeg*Nz; i++)HconnT[i] = -1;
	for (i = 0; i < N1 - K1; i+=z)
	{
		for (j = 0; j < ncheck; j++)
		{
			if (Hconn[i*ncheck + j] < 0)break;
			pcol = HconnT + (Hconn[i*ncheck + j])/z * maxdeg;
			while (*pcol >= 0)pcol++;
			*pcol = i/z;
		}
	}

	printf("MakeEncodable: Permuting PCM to lower triangular form...\n");
	// Contruct row and column permutations minimizing gap
	int *pr = new int[Nz - Kz];
	int* prinv = new int[Nz - Kz];
	for (i = 0; i < Nz - Kz; i++)prinv[i] = pr[i] = i;

	int *pc = new int[Nz];
	int* pcinv = new int[Nz];
again:
	for (i = 0; i < Nz; i++)pc[i] = pcinv[i] = i;// Contain the index of the i-th column in the permuted matrix


	gap = 0;
	for (n = 0; n < Nz - Kz - gap; n++)
	{
		for (i = 0; i < Nz - Kz; i++)pr[prinv[i]] = i;

		// (2) find column with smallest number of '1' above g 
		// prinv[i] is the index in the new matrix of original row i
		// pr[i] is the index in the original matrix of new row i
		ncmin = 10000;
		degmin = maxdeg;
		for (n1 = n; n1 <Nz; n1++)// check all available columns
		{
			pcol = HconnT + pc[n1] * maxdeg; // Point to the column pc[n1] is the index in the original matrix of n1 column in the new
			nc = 0;
			for (k = 0; k < maxdeg; k++)// index of rows
			{
				if (pcol[k] < 0)break; // Row number
				if (prinv[pcol[k]] >= n && prinv[pcol[k]]< Nz - Kz - gap)nc++;
			}
			//			printf("%d\t%d\n", n1, nc);
			if (nc != 0)
			{
				if (nc < ncmin)
				{
					degmin = k;
					ncmin  = nc;
					cbest  = n1;
				}
				else if (nc == ncmin && k<degmin )
				{
					cbest = n1;
					degmin = k;
				}
			}
		}

		// Swap column N1-n1 with cbest
//		printf("n=%d,  lowest degree at %d (%d)\n", n, cbest, ncmin);
		if (n != cbest)
		{
			temp = pc[n];
			pc[n] = pc[cbest];
			pc[cbest] = temp;
			pcinv[pc[n]] = n;
			pcinv[pc[cbest]] = cbest;
			printf("%3d: Swap Col %3d <-> %3d \n",n,pc[n], pc[cbest]);
		}

		for (i = 0; i < Nz; i++)pcinv[pc[i]] = i;


		pcol = HconnT + pc[n] * maxdeg; // Point to the active column
		for (k = l = 0; k < maxdeg; k++)
		{
			if (pcol[k] < 0)break;
			if (prinv[pcol[k]] >= n && prinv[pcol[k]]< Nz - Kz - gap)
			{
				if (l == 0)
				{
					temp = prinv[pcol[k]];
					prinv[pcol[k]] = n;
					prinv[pr[n]] = temp;
				}
				else
				{
					temp = prinv[pcol[k]];
					prinv[pcol[k]] = Nz - Kz - 1 - gap;
					prinv[pr[Nz - Kz - 1 - gap]] = temp;
					gap++;
				}
				// printf("Row %d is swapped with %d\n", prinv[pr[g + l]],prinv[pcol[k]]);
				l++;
			}
		}
	}

	for (i = 0; i < Nz - Kz; i++)pr[prinv[i]] = i;
	printf("MakeEncodable: Residual Gap = %d x %d.\n", gap,z);
//	goto again;


	///////////// Print Lower triangular matrix ///////////////////////
	//int* HH = (int*)calloc((Nz - Kz)*Nz, 4);
	//int row;
	//for (i = 0; i < Nz; i++)
	//{
	//	for (j = 0; j < maxdeg; j++)
	//	{
	//		row = HconnT[i * maxdeg + j];
	//		if (row < 0)break;// Contain row number (index of check nodes)
	//		if(prinv[row]<Nz-Kz-gap)HH[(Nz-Kz-1-gap-prinv[row]) * Nz + (Nz-1-pcinv[i])] = 1;
	//		else                    HH[(prinv[row])             * Nz + (Nz-1-pcinv[i])] = 1;
	//	}
	//}
	//int deg;
	//for (j = 0; j < Nz; j++)
	//{
	//	if (j == K1)printf("|");
	//	printf("%2d", nconnvar[pc[Nz-1-j]]);
	//}
	//printf("\n");
	//for (i = 0; i < Nz - Kz; i++)
	//{
	//	deg = 0;
	//	for (j = 0; j<Nz; j++)
	//	{
	//		if (j == K1)printf("|");
	//		printf("%2d", HH[i*Nz + j]);
	//		deg += HH[i*Nz + j];
	//	}
	//	printf("\t%d\n",deg);
	//}
	//delete[] HH;
	////////////////////////////////////
	delete[] HconnT;
	HconnT = (int*)calloc(ncheck*(N1 - K1),sizeof(int));
	int*pH;
	int l1, l2, col, shift, jmax;
	if (pattern != 0)// Update puncturing pattern
	{
		int* newpat = new int[N1];
		for (col = 0; col < Nz; col++)
		{
			for (k = 0; k < z; k++)
			{
				newpat[(Nz - 1 - pcinv[col])*z+k] = pattern[col*z + k];
			}
		}
		delete[] pattern;
		pattern = newpat;
	}
	for (i = 0; i < Nz - Kz; i++)
	{
		//		pH = HconnT + (N1 - K1 - 1 - prinv[i]) * ncheck;
		if (prinv[i] < Nz - Kz - gap)
		{
			pH = HconnT + (Nz - Kz  - gap - prinv[i]- 1) * ncheck*z;
		}
		else
		{
			pH = HconnT + (prinv[i]) * ncheck*z;
		}

		for (j = 0; j < ncheck; j++)
		{
			col = Hconn[i*z*ncheck + j];
			if (col < 0)break;
			shift = (col % z);
			col = (Nz - 1 - pcinv[col / z]); // Flip matrix 
			for (k = 0; k < z; k++)
			{
				pH[k*ncheck + j] = col*z + ((shift + k) % z);
			}
		}
		if(j<ncheck)for (k = 0; k < z; k++)pH[k*ncheck + j] = -1;
		// Sorting matrix and place identity shift at the right
		shift = 0;
		jmax = j;
		for (k = 0; k < z; k++)
		{

			for (l1 = 0; l1 < jmax - 1; l1++)
			{
				for (l2 = l1 + 1; l2 < jmax; l2++)
				{
					if (pH[l1] > pH[l2]) { temp = pH[l1]; pH[l1] = pH[l2]; pH[l2] = temp; }
				}
			}
			if (k == 0)shift = (pH[jmax-1]%z);
			if(shift!=0)
				for (l1 = 0; l1 < jmax; l1++)pH[l1] = (pH[l1] / z)*z + (pH[l1]+ z - shift) % z;
			for (j=jmax; j < ncheck; j++)pH[j] = -1;
			pH += ncheck;
		}
	}

	delete[] pc;
	delete[] pcinv;
	delete[] pr;
	delete[] prinv;
	delete[] Hconn;
	Hconn = HconnT;
	delete[] nconnvar;
	delete[] gg; gg = 0;
	if (gap > 0)
	{
		printf("MakeEncodable: Computing Gap matrix...");
		ComputeGapMatrix(gap,z);
		// Sort matrix again
		for (i = 0; i < Nz - Kz; i++)
		{
			pH = Hconn + i * ncheck*z;

			for (k = 0; k < z; k++)
			{
				for (l1 = 0; l1 < ncheck - 1; l1++)
				{
					if (pH[l1] < 0)break;
					for (l2 = l1 + 1; l2 < ncheck; l2++)
					{
						if (pH[l2] < 0)break;
						if (pH[l1] > pH[l2]) { temp = pH[l1]; pH[l1] = pH[l2]; pH[l2] = temp; }
					}
				}
				pH += ncheck;
			}
		}
		printf("Completed.\n");
		gap *= z;
		SetEncodingMethod(5);
	}
	else
	{
		SetEncodingMethod(0);
	}
}
void LDPC_Encoder::ComputeGapMatrix(int g,const int z)
{
	// The Hconn matrix should be in the best possible lower triangular form, with residual gap g
	// The entries of each row in Hconn must be sorted
	// Computes the matrix to be used for the computation of gap bits

	g *= z;
	int i, j, ii,k;
	int m = N1 - K1;
	int mg = m - g;

	int* AAA = (int*)calloc((N1 - K1)*K1, sizeof(int));
	int* BBB = (int*)calloc((N1 - K1)*g, sizeof(int));
	int* TTT = (int*)calloc(mg*mg, sizeof(int));
	int *EEE = (int*)calloc(g*mg, sizeof(int));
		
	printf("Computes matrices A|C), B|D), T \n");
	//	Computes matrices A|C), B|D), T and E
	for (i = 0; i < N1-K1; i++)
	{
		if (i < mg)
		{
			for (j = 0; j < ncheck; j++)
			{
				if (Hconn[i*ncheck + j]>=0)
				{
					if (Hconn[i*ncheck + j] < K1)			AAA[i*K1+Hconn[i*ncheck + j]]          = 1;
					else if (Hconn[i*ncheck + j] < K1 + g)	BBB[i*g+Hconn[i*ncheck + j] - K1]      = 1;
					else									TTT[i*mg+Hconn[i*ncheck + j] - K1 - g] = 1;
				}
				else break;
			}
		}
		else
		{
			for (j = 0; j < ncheck; j++)
			{
				if (Hconn[i*ncheck + j] >= 0)
				{
					if (Hconn[i*ncheck + j] < K1)			AAA[i*K1 + Hconn[i*ncheck + j]] = 1;  //C
					else if (Hconn[i*ncheck + j] < K1 + g)	BBB[i*g + Hconn[i*ncheck + j] - K1] = 1; //D
					else									EEE[(i - mg)*mg + Hconn[i*ncheck + j] - K1 - g] = 1;
				}
				else break;
			}

		}
	}	

	// T-1
	printf("Inverse of lower triangular %d....", mg);
	int* INVTTT = (int*)calloc(mg*mg, sizeof(int));
	int* p1,*p2;
	for (i = 0; i < mg; i++)// Inverse of lower triangular
	{
		for (j = 0; j < i; j++)
		{
			if (TTT[i*mg + j] != 0)
			{
				p1 = INVTTT + i*mg;
				p2 = INVTTT + j*mg;
				for (k = 0; k <= j; k++)p1[k] ^= p2[ k];
			}
		}
		INVTTT[i*mg + i] = 1;
	}
//	int* check = mmult(TTT, INVTTT,mg,mg,mg);
	free(TTT);
	printf("Done\n");


	// (ET-1|I) is a matrix gx(mg,g) the right part is the identity
	int* PRE = (int*)calloc(g*m, sizeof(int));
	int t;
	for (i = 0; i < g; i++)PRE[i*m + mg + i] = 1;
	printf("MMult %dx%dx%d\n", g, mg, mg);
	for (i = 0; i < g; i++)
	{
		p1 = EEE    + i*mg;
		for (ii = 0; ii < mg; ii++)
		{
			if (p1[ii] == 0)continue;
			p2 = INVTTT + ii*mg;
			for (j = 0; j < mg; j++)PRE[i*m + j] ^= p2[j];
		}
	}

	// Compute matrices (HHH =ET-1A+C and PHI=(-ET-1B+D), and store them in H and PHI)
	printf("MMult %dx%dx%d\n", g,m,K1);
	int* HHH = mmult(PRE, AAA, g, m, K1);
	printf("MMult %dx%dx%d\n", g,m,g);
	int* PHI = mmult(PRE, BBB, g, m, g); 
	free(EEE);
	free(PRE);
	free(AAA);
	free(BBB);

	int* INVPHI = 0;
	int col1, col2, temp,s;
	int seed = 2347892;
	int jj;
	int *perm = new int[(K1 + g)/z];
	int *perminv = new int[(K1 + g)/z];
	for (i = 0; i < (K1 + g)/z; i++)perm[i] = i;
	// Randomly swap macro-columns until we find an invertible matrix  
	for (j = 0; j < 100; j++)
	{
		printf("Inverse of dimension %d....", g);
		INVPHI = Inverse_bin(PHI, g);
		printf("Done\n.");
		if (INVPHI != 0)break;
		col1 = unif_int(seed) % (K1/z);
		col2 = unif_int(seed) % (g/z);
		// Swap macro-columns
		for (ii = 0; ii < g; ii++)
		{
			for (jj = 0; jj < z; jj++)
			{
				temp                   = HHH[ii*K1 + col1*z + jj];
				HHH[ii*K1 + col1*z + jj]= PHI[ii*g + col2*z + jj];
				PHI[ii*g + col2*z + jj]= temp;
			}

		}
		temp = perm[col1];
		perm[col1] = perm[(K1/z)+col2];
		perm[(K1/z)+col2] = temp;
	}
	//printf("------ PHI -------\n");
	//for (i = 0; i < g*g; i++)
	//{
	//	printf("%1d", PHI[i]);
	//	if (i%g == g - 1)printf("\n");
	//}


	if (INVPHI != 0)
	{
		/* Optional Check */
		int ttt = 0;
		for (i = 0; i < g; i++)
		{
			for (j = 0; j < g; j++)
			{
				ttt = 0;
				for (k = 0; k < g; k++)
				{
					ttt ^= PHI[i*g + k] * INVPHI[k*g + j];
				}
				if (i == j && ttt != 1)printf("ComputeGapMatrix: Error inversion.\n");
				if (i != j && ttt != 0)printf("ComputeGapMatrix: Error inversion.\n");
			}
			//		printf("\n", ttt);
		}

		for (i = 0; i < (K1 + g)/z; i++)perminv[perm[i]] = i;
		//printf("------ INV--------\n");
		//for (i = 0; i < g*g; i++)
		//{
		//	printf("%1d", INVPHI[i]);
		//	if (i%g == g - 1)printf("\n");
		//}

		if (pattern != 0)// Update Puncturing pattern
		{
			for (i = 0; i < (K1 + g)/z; i++)
			{
				if (perminv[i] == i)continue;
				for (k = 0; k < z; k++)
				{
					temp                      = pattern[i*z+k]; 
					pattern[i*z + k]          = pattern[perminv[i] *z + k];
					pattern[perminv[i]*z + k] = temp;
				}
			}
		}
		for (i = 0; i < (N1-K1)/z; i++)// row index
		{
			for (j = 0; j < ncheck; j++)
			{
				temp = Hconn[i*z*ncheck + j];
				if (temp < 0)break;
				s    = temp % z;
				temp /= z;  //macrocolumn
				if (temp >= (K1 + g) / z)continue;
				if (perminv[temp] == temp)continue;
				for (k = 0; k < z; k++)
				{
					Hconn[(i*z+k)*ncheck  + j] = perminv[temp]*z + (s+k)%z;
				}
			}
		}
		gg = mmult(INVPHI, HHH, g, g, K1);
		free(INVPHI);
	}
	else
	{
		gg = 0;
	}
	delete[]perm;
	delete[]perminv;
	free(HHH);
	free(PHI);
}
void LDPC_Encoder::SetRegularRandom2(const int K,const int N,double avdeg)
{
	int i,j,jj,ipick,seed,temp;
	this->K=K;
	this->N=N1=N;
	
	seed=712991519;

	if(avdeg>0)this->ncheck=(int)((double)(avdeg*N)/(N-K)+0.9999);
	else	   this->ncheck=(int)((double)(2*N)/(N-K)+0.9999);
	delete[] Hconn;
	Hconn=new int[this->ncheck*(N-K)];

	int *vect=new int[N];
	int l;
	bool froze;
	l=0;
	for(i=0;i<N;i++)vect[i]=i;
	for(i=0;i<N-K;i++)
	{
		froze=false;
		for(j=0;j<ncheck;j++)
		{
			if (j < ncheck - 1)
			{
				do {

					if(froze)ipick=l;
					else   	 ipick=l + (unif_int(seed)%(N-l));
				} while (vect[ipick] >= K + i);

			}
			else
			{
				ipick = 0;
				while (vect[ipick] != K +i)ipick++;
			}

			temp=vect[l];
			vect[l]=vect[ipick];
			vect[ipick]=temp;
			Hconn[i*ncheck+j]=vect[l];
			l++;
			if(l==N)
			{
				l=0;
				froze=true;
			}
		}
		// Sort [optional]
		for(j=0;j<ncheck-1;j++)
		{
			for(jj=j+1;jj<ncheck;jj++)
			{
				if(Hconn[i*ncheck+j]>Hconn[i*ncheck+jj])
				{
					temp=Hconn[i*ncheck+jj];
					Hconn[i*ncheck+jj]=Hconn[i*ncheck+j];
					Hconn[i*ncheck+j]=temp;
				}
			}
		}
	}
	delete[] vect;
}
void LDPC_Encoder::SetRegularRandom(const int K,const int N,double avdeg)
{
	int i,j,jj,ipick,seed,temp;
	this->K=K1=K;
	this->N=N1=N;
	
	seed=712991519;

	if(avdeg>0)this->ncheck=(int)((double)(avdeg*N)/(N-K)+0.9999);
	else	   this->ncheck=(int)((double)(2*N)/(N-K)+0.9999);
	delete[] Hconn;
	Hconn=new int[this->ncheck*(N-K)];

	int *vect=new int[N];
	int l;
	bool froze;
	l=0;
	for(i=0;i<N;i++)vect[i]=i;
	for(i=0;i<N-K;i++)
	{
		froze=false;
		for(j=0;j<ncheck-2;j++)
		{
			if(froze)ipick=l;
			else   	 ipick=l + (unif_int(seed)%(K-l));

			temp=vect[l];
			// Check if it has been picked by previous check (avoid 
			if(i>0)
			{
				for(jj=0;jj<ncheck;jj++)
				{
					if(Hconn[(i-1)*ncheck+jj]==temp)
					{
						printf("warning.\n");
					}
				}
			}
			vect[l]=vect[ipick];
			vect[ipick]=temp;
			Hconn[i*ncheck+j]=vect[l];
			l++;
			if(l==K)
			{
				l=0;
				froze=true;
			}
		}
		Hconn[i*ncheck+ncheck-2] = K+i-1;
		Hconn[i*ncheck+ncheck-1] = K+i;
		// Sort [optional]
		for(j=0;j<ncheck-2;j++)
		{
			for(jj=j+1;jj<ncheck-1;jj++)
			{
				if(Hconn[i*ncheck+j]>Hconn[i*ncheck+jj])
				{
					temp=Hconn[i*ncheck+jj];
					Hconn[i*ncheck+jj]=Hconn[i*ncheck+j];
					Hconn[i*ncheck+j]=temp;
				}
			}
		}
	}
	delete[] vect;
}
void LDPC_Encoder::SetIrregularRandom(const int K,const int N,const int ndv,const double* pdv,bool edgeper)
{
	int i,j,k,ipick,seed,temp;
	this->K=K;
	this->N=N1=N;
	int nones;
	double tot;

	int* ldegn=new int[ndv]; // Number of nodes of degree i

	if(edgeper)
		for(i=0,tot=0;i<ndv;i++){tot+=pdv[2*i+1]/pdv[2*i];}
	else
		for(i=0,tot=0;i<ndv;i++){tot+=pdv[2*i+1];}


	nones=0;
	int m=0;
	for(i=0;i<ndv-1;i++)
	{
		if(edgeper)	m+=ldegn[i]= (int)(N*pdv[2*i+1]/pdv[2*i]/tot);
		else		m+=ldegn[i]= (int)(N*pdv[2*i+1]/tot);
		nones+=ldegn[i]*(int)pdv[2*i];
	}
	ldegn[ndv-1]=N-m;
	nones+=ldegn[ndv-1]*(int)pdv[2*(ndv-1)];
	int ncheck=(nones/(N-K))+1;

	// Initialize vectors
	int *vect=new int[nones];
	int l,l2;
	l=l2=0;
	for(i=0;i<ndv;i++)
	{
		for(j=0;j<ldegn[i];j++)
		{
			for(k=0;k<pdv[2*i];k++)
			{
				vect[l]=l2;
				l++;
			}
			l2++;
		}
	}

	// Random permutation
	seed=712991519;
	for(l=0;l<nones;l++)
	{
		ipick=l + (unif_int(seed)%(nones-l));
		temp=vect[l];
		vect[l]=vect[ipick];
		vect[ipick]=temp;
	}

	delete[] Hconn;
	Hconn=new int[(N-K)*ncheck];

	l=0;
	for(j=0;j<ncheck;j++)
	for(i=0;i<(N-K);i++)
	{
		if(l<nones)Hconn[i*ncheck+j]=vect[l++];
		else	   Hconn[i*ncheck+j]=-1;
	}
	delete[] ldegn;
	delete[] vect;

	this->ncheck=ncheck;
}
void LDPC_Encoder::SetRightIrregularRandom(const int K,const int N,const double avdeg,const int ndc,const int* dc,const double* pdc,bool edgeper)
{
	int i,j,k,ipick,seed,temp;
	this->K=K;
	this->N=N1=N;
	int nones;
	double tot,tot2;

	int* rdegn=new int[ndc]; // Number of nodes of degree i
	int* dct=new int[ndc];

	double avedc=avdeg/(1.-(double)K/N);


//	double alpha=avedc*tot/tot2;
//	double avedc=avdeg/(1.-(double)K/N);
	//for(i=0,tot=tot2=0;i<ndc;i++){tot+=pdc[i];tot2+=pdc[i]*dc[i];}

	//int ncheck=0;
	//int m=0;
	//nones=0;
	//for(i=0;i<ndc;i++)
	//{
	//	dct[i]=(int)(alpha*dc[i]+0.5);
	//	if(i<ndc-1)rdegn[i]=(int)(pdc[i]/tot*(N-K)+0.5); // Number of check node with degree di;
	//	else rdegn[ndc-1]=N-K-m;
	//	m+=rdegn[i];
	//	nones+=rdegn[i]*dct[i];
	//	if(dct[i]>ncheck)ncheck=dct[i];
	//}

	for(i=0,tot=tot2=0;i<ndc;i++){tot+=pdc[i];tot2+=pdc[i]*(dc[i]-dc[0]);}
	double alpha=(avedc-dc[0])*tot/tot2;
	int ncheck=0;
	int m=0;
	nones=0;
	for(i=0;i<ndc;i++)
	{
		dct[i]=dc[0]+(int)(alpha*(dc[i]-dc[0])+0.5);
		if(dct[i]>ncheck)ncheck=dct[i];
		if(i<ndc-1) rdegn[i]=(int)(pdc[i]/tot*(N-K)+0.5); // Number of check node with degree di;
		else		rdegn[i]=N-K-m;
		m		+=rdegn[i];
		nones	+=rdegn[i]*dct[i];
	}
	if(nones<2*N)
	{
		printf("SetRightIrregularRandom: increase the average check degree above %f\n",(2.*N)/(N-K));
		exit(1);
	}
	for(i=0;i<ndc;i++)
	{
		printf("%d\t%d\t%d\n",i,dct[i],rdegn[i]);
	}


	int ldeg  = (nones/N)+1;
	printf("Irregular LDPC density: L=%f  R=%f\n",(double)nones/N, (double)nones/(N-K));


	// Initialize vectors
	int *vect=new int[nones];
	int l,l2;
	for(l=0;l<nones;l++)vect[l]=l%N;

	// Random permutation
	seed=712991519;
	for(l=0;l<nones;l++)
	{
		ipick=l + (unif_int(seed)%(nones-l));
		temp=vect[l];
		vect[l]=vect[ipick];
		vect[ipick]=temp;
	}

	delete[] Hconn;
	Hconn=new int[(N-K)*ncheck];

	l2=l=0;
	for(i=0;i<ndc;i++)
	{
		for(j=0;j<rdegn[i];j++)
		{
			for(k=0;k<dct[i];k++)Hconn[l*ncheck+k]=vect[l2++];
			for(;k<ncheck;k++)	   Hconn[l*ncheck+k]=-1;
			l++;
		}
	}
	delete[] dct;
	delete[] rdegn;
	delete[] vect;

	this->ncheck=ncheck;
}
//void LDPC_Encoder::SetLeftIrregularRandom(const int K, const int N, 
//	const int ndv, //! Number of terms
//	const int* dv,  // degrees
//	const double* pdv,	// probabilities
//	bool edgeper, bool PEG)
//{
//	int i, j, k, ipick, seed, temp;
//	this->K = K1 = K;
//	this->N = N1 = N;
//	int nones;
//	double tot;
//
//	int* ldegn = new int[ndv]; // Number of nodes of degree dv[i] (1 is of degree 1)
//	int* ldeg  = new int[N]; 
//	int* rdeg  = new int[N-K]; // Number of nodes of degree i (1 is of degree 1)
//
//	if (edgeper)
//		for (i = 0, tot = 0; i<ndv; i++) { tot += pdv[i] / dv[i]; }
//	else
//		for (i = 0, tot = 0; i<ndv; i++) { tot += pdv[i]; }
//
//
//	nones = 0;
//	int m = 0;
//	for (i = 0; i<ndv - 1; i++)
//	{
//		if (edgeper)ldegn[i] = (int)(N*pdv[i] / dv[i] / tot);
//		else        ldegn[i] = (int)(N*pdv[i] / tot);
//		for (j = 0; j < ldegn[i]; j++)ldeg[m + j] = dv[i];
//		m += ldegn[i];
//		nones += ldegn[i] * (int)dv[i];
//	}
//	ldegn[ndv - 1] = N - m;
//	for (j = 0; j < ldegn[ndv - 1]; j++)ldeg[m + j] = dv[ndv-1];
//	nones += ldegn[ndv - 1] * dv[ndv - 1];
//
//	for (i = 0; i<ndv; i++)printf("%d\t%d (%f)\n", dv[i], ldegn[i], pdv[i]);
//
//
//	// Compute number of check with degree rmin and rmin+1
//	printf("%d\t%f\n", nones, (double)nones / N);
//	double etar = (double)nones / (N - K);
//	int    rmin = (int)etar;
//	int R[4];
//	R[1] = nones % (N - K);  //c
//	R[0] = (N - K) - R[1];
//	for (i = 0; i < N - K; i++)
//	{
//		if(i<R[0])rdeg[i] = rmin;
//		else	  rdeg[i] = rmin+1;
//	}
//
//
//	ncheck = rmin + 1;
//
//	if (PEG)// Matrix construction with PEG
//	{
//		LDPC_Encoder::PEG(K, N, nones, ldeg,rdeg);
//	}
//	else
//	{
//
//	}
//	delete[] ldeg;
//	delete[] rdeg;
//
//	this->ncheck = ncheck;
//}
void LDPC_Encoder::SetLeftIrregularRandom(
	const int K, 
	const int N, 
	const int ndv, 
	int* dv, 
	const double* pdv, 
	const int Z, 
	const bool edgeper,
	const bool PEG)
{
	int i, j, k, ipick, seed, temp;
	this->K = K1 = K*Z;
	this->N = N1 = N*Z;
	this->z = Z;
	int nones;
	int l, l2;

	double tot;

	int* ldegn = new int[ndv]; // Number of nodes of degree i (1 is of degree 1)

	if (edgeper)
		for (i = 0, tot = 0; i<ndv; i++) { tot += pdv[i] / dv[i]; }
	else
		for (i = 0, tot = 0; i<ndv; i++) { tot += pdv[i]; }


	int max,T;
	bool left = true;
	if (left)
	{
		max = (N - K);
		T = N;
	}
	else
	{
		max = N;
		T = N-K;
	}

	nones = 0;
	int m = 0;
	double cum = 0.;
	int d;
	for (i = 0; i < ndv; i++)
	{
		if (dv[i] > max)
		{
			printf("Clipping degree %d to max (%d)\n", dv[i], N - K);
			dv[i] = max;
		}
	}
	for (i = 0; i<ndv - 1; i++)
	{
		cum += pdv[i] / tot;
		d = (int)(cum*T + 0.5) - m;
		ldegn[i] = (int)(T*pdv[i] / dv[i] / tot+0.5);		
		if (edgeper)	m += ldegn[i];
		else
		{
			m += ldegn[i] = d;
			nones += ldegn[i] * (int)dv[i];
		}
	}
	ldegn[ndv-1] = T - m;
	nones += ldegn[ndv-1 ] * dv[ndv - 1];

	printf("Finite size (%dx%d) degree  distribution\n",N,Z);
	for (i = 0; i<ndv; i++)
		if(ldegn[i]>0)printf("N[%d]=\t%d \t(%f)\n", dv[i], ldegn[i],(double)ldegn[i]/T);


	printf("Edges=%d\teta=%f\n", nones*Z, (double)nones / N);


	if(PEG)// Matrix construction with PEG
	{
		int sglConcent = 1;  // Set to 0 for forcing concentration 1 otherwise
		int targetGirth = 100;
		bool verbose = true;
		BigGirth *bigGirth;
		int* degseq = new int[T];
		for (i =l= 0; i < ndv; i++)
		{
			for (j = 0; j < ldegn[i]; j++) degseq[l++] = dv[i];
		}
		printf("Running PEG optimization (BigGirth...");
		bigGirth = 0;

		bigGirth=new BigGirth((N - K)*Z, N*Z, degseq, sglConcent, targetGirth, Z,!left);

		//bigGirth->writeToFile_Hmatrix("H.txt");
	//	bigGirth->writeToFile_Hcompressed("Hcomp.txt",Z);
		printf("Completed.\n");
		delete[] Hconn; Hconn = 0;
//		delete[] Gconn; Gconn = 0;
		ncheck = 0;
		Hconn=bigGirth->Build_Hcompressed(ncheck);
		delete bigGirth;
		delete[] degseq;
		delete[] ldegn;
		return;
	}
	else  // Random Construction
	{
		int *vect = new int[nones];
		int *ulim = new int[(N - K) + 1];
		l = l2 = 0;
		for (i = ndv - 1; i >= 0; i--)
		{
			for (j = 0; j < ldegn[i]; j++)
			{
				for (k = 0; k < dv[i]; k++)
				{
					vect[l] = l2;
					l++;
				}
				l2++;
				if (l2 >= K)ulim[l2 - K] = l;
			}
		}
		nones--;
		double etar = (double)nones / (N - K);
		int    rmin = (int)etar;
		int R[4];
		R[1] = nones % (N - K);  //c
		R[0] = (N - K) - R[1];

		ncheck = rmin + 1;

		int deg;

		delete[] Hconn;
		Hconn = new int[(N - K)*ncheck];

		// Random permutation constrained to avoid parallel edges
		seed = 712991519;
		int row, col;
		row = col = 0;
		for (row = l = 0; row < N - K; row++)
		{
			if (row < R[0])deg = rmin;
			else			deg = rmin + 1;

			for (col = 0; col < deg - 1; col++)
			{
				if ((ulim[row] - l)>0)ipick = l + (unif_int(seed) % (ulim[row] - l));
				else			   ipick = l;
				temp = vect[l];
				vect[l] = vect[ipick];
				vect[ipick] = temp;
				Hconn[row*ncheck + col] = vect[l];
				l++;
			}
			ipick = ulim[row];
			temp = vect[l];
			vect[l] = vect[ipick];
			vect[ipick] = temp;
			Hconn[row*ncheck + col] = vect[l];  // Should be K+row
			if (vect[l] != K + row)
			{
				printf("Warning");
			}
			//		Hconn[row*ncheck + deg-1] = K+row;
			if (col < rmin)Hconn[row*ncheck + col + 1] = -1;
			l++;
		}

		delete[] ldegn;
		delete[] ulim;
		delete[] vect;
		this->ncheck = ncheck;
	}



}



void LDPC_Encoder::RunM(const int ntics,const int* inp, int* out)
{
	int i,j,ii,temp;
	int* h,*g;
	for(ii=0;ii<ntics;ii++)
	{

		switch(enctype)
		{
		case 0: // Encoding with parity check in lower triangular form with back substitution
			h=Hconn;
			for(i=0;i<K;i++)out[i]=inp[i];// Systematic bits
			for(;i<N;i++)
			{
				temp=0;
				for(j=1;j<ncheck;j++)
				{
					if(h[j]==-1)break;	
					temp += out[h[j-1]];
				}
				temp%=M;
				if(temp>0)temp=M-temp;
				out[h[j-1]]=temp;
				h+=ncheck;
			}
			break;

		case 1:// Encoding with Generator matrix
			int l,bit;
			g=G;
			for(i=0;i<K;i++)out[i]=inp[i];// Systematic bits
			for(j=l=0;j<N-K;j++)
			{
				for(i=bit=0;i<K;i++)
				{
					bit += g[i]*inp[i];
				}
				bit%=M;
				if(bit!=0)bit=M-bit;
				out[K+j]=bit;
				g+=K;
			}
			break;

		default: 
			printf("Encoding type %d not available for mod-M encoding\n",enctype);
			exit(1);
			break;
		}
		out+=N;
		inp+=K;
	}
}
void LDPC_Encoder::Run(const int ntics,const int* inp, int* outt)
{
	int i,j,ii,l,temp;
	int* h;
	int* out;
	int* ggg = gg;
	K1=K+shorten;
	for(ii=0;ii<ntics;ii++)
	{
		if(period==0 && shorten==0)// On place encoding
		{
				out = outt;
				for(i=0;i<K;i++)out[i]=inp[i];// Systematic bits
		}
		else	// Use buffer
		{
			out = buff;
			if (lastsh == 0)
			{
				for(i=0;i<shorten;i++)	out[i]=0;// shortened bits
				for(i=0;i<K;i++)		out[i+shorten]=inp[i];// Systematic bits

			}
			else
			{
				for(i=0;i<K;i++)		out[i]=inp[i];// Systematic bits
				for(i=0;i<shorten;i++)	out[K+i]=0;// shortened bits
			}
		}

		switch(enctype)
		{
		case 0: // Default Encoding: Encoding with parity check in lower triangular form with back substitution
			h=Hconn;
			for(i=K1;i<N1;i++)
			{
				temp=0;
				for(j=1;j<ncheck;j++)
				{
					if(h[j]==-1)break;	
					temp ^= out[h[j-1]];
				}
				out[h[j-1]]=temp;
				h+=ncheck;
			}
			break;



		case 1:// Encoding with Generator matrix
			int l, bit;
			h = G;
			for (j = l = 0; j<N1 - K1; j++)
			{
				for (i = bit = 0; i<K1; i++)
				{
					bit ^= h[i] & out[i];
				}
				out[K1 + j] = bit;
				h += K1;
			}
			break;
		case 2: // Method 1a in 802.16e-2005  Matrix in H in circular block form

				// Initial vector v(0)
			for (i = K1; i<K1 + z; i++)out[i] = 0;
			// Initial vector v(0)
			for (j = 0; j<N1 - K1; j++)
			{
				h = Hconn + ncheck*j;
				l = 0;
				while (h[l]<K1 && h[l] != -1)
				{
					out[K1 + ((j + x) % z)] ^= out[h[l]]; // Initializing bits
					l++;
				}
			}

			for (; i<N1; i++)// z is the expanding factor
			{
				out[i] = 0;
				if (i >= K1 + 2 * z)
				{
					out[i] ^= out[i - z];
				}
				h = Hconn + ncheck*(i - K1 - z);
				l = 0;
				while (h[l]< K1 + z && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
			}
			break;

		case 3:// Encoding with Compressed Generator matrix
			h = Gconn;
			for (j = l = 0; j<N1 - K1; j++)
			{
				for (i = bit = 0; i<nmaxvar; i++)
				{
					if (h[i] < 0)break;
					bit ^= out[h[i]];
				}
				out[K1 + j] = bit;
				h += nmaxvar;
			}
			break;

		case 4: // Method NG-LTE Qualcomm
			// Initial vector v(0) sum of all parity checks
			for(i=K1+z;i<K1+2*z;i++)out[i]=0;
			for(j=0;j<(ccore-1)*z;j++)
			{
				h=Hconn + ncheck*j;
				l=0;
				while(h[l]<K1 && h[l] != -1)
				{
					out[K1 + z+ (j%z)] ^= out[h[l]]; // Initializing bits
					l++;
				}
			}
			if (aaa!=0 || hhh!= 0)// Special encoding for family 1
			{

				int* zz2 = out + K1;
				int* zz1 = out + K1 + z;
				if (hhh > 0)
				{
					for (i = 0; i < z; i++)// z2=z1*Q(x)=1+x^5h+x^12h
					{
						//zz2[i] = zz1[i]
						//	^ zz1[(1000 * z + i - 5 * hhh) % z]
						//	^ zz1[(1000 * z + i - 12 * hhh) % z];
						zz2[i] = zz1[i]
							^ zz1[(i + 5 * hhh) % z]
							^ zz1[(i + 12 * hhh) % z];
					}
					for (i = 0; i < z; i++)// Q(x*x) 1+x^10h+x^24h x^-a)
					{
						//zz1[i] = zz2[(1000 * z + i + aaa) % z]
						//	^ zz2[(1000 * z + i - 10 * hhh + aaa) % z]
						//	^ zz2[(1000 * z + i - 24 * hhh + aaa) % z];
						zz1[i] = zz2[(i       + z - aaa) % z]
							^ zz2[(i + 10 * hhh + z - aaa) % z]
							^ zz2[(i + 24 * hhh + z - aaa) % z];
					}
				}
				else
				{
					for (i = 0; i < z; i++)zz2[i] = zz1[i];
					for (i = 0; i < z; i++)zz1[i] = zz2[(z-aaa) % z];
				}
			}
			// Core
			for (i=K1 + 2 * z; i< K1 + ccore*z; i++)// Core
			{
				out[i] = 0;
				h = Hconn + ncheck*(i - K1 - 2*z);
				l = 0;
				while (h[l]< i && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
				if (h[l] != i)
				{
					printf("Error\n");
				}
			}

			// Parity check of first two columns
			for (i = K1; i < K1 + z; i++)
			{
				out[i] = 0;
				h = Hconn + ncheck*((ccore-1)*z+i-K1);
				l = 0;
				while (h[l] < i && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
				if (h[l] != i)
				{
					printf("Error\n");
				}
			}

			// HARQ section
			for (i=K1 + ccore*z; i< N1; i++)// Core
			{
				out[i] = 0;
				h = Hconn + ncheck*(i - K1);
				l = 0;
				while (h[l]< i && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
				if (h[l] != i)
				{
					printf("Error\n");
				}
			}
			break;

		case 5: //Encoding with Richardson Algorithm
			if (gg == 0)exit(1);
			h = gg;
			// Generate gap bits
			for (i = K1; i < K1 + gap; i++)
			{
				out[i] = 0;
				for (j = 0; j < K1; j++)
				{
					if (h[j] != 0)out[i] ^= out[j];

				}
				h += K1;
			}
			h = Hconn;
			for (i = K1 + gap; i<N1; i++)
			{
				temp = 0;
				for (j = 1; j<ncheck; j++)
				{
					if (h[j] <0)break;
					if (h[j - 1] > i)printf("Error.\n");
					temp ^= out[h[j - 1]];
				}
				if (h[j - 1] != i)
				{
					printf("Error.\n");
				}

				out[i] = temp;
				h += ncheck;
			}
			break;
		case 6:
			// Initial vector v(0) sum of all parity checks	
			for (i = K1; i<K1 + z; i++)out[i] = 0;
			for (j = 0; j<ccore*z; j++)
			{
				h = Hconn + ncheck*j;
				l = 0;
				while (h[l]<K1 && h[l] != -1)
				{
					out[K1 + ((j+aaa)%z)] ^= out[h[l]]; // Initializing bits
					l++;
				}
			}
			// Core
			for (i = K1 + z; i< K1 + ccore*z; i++)// Core
			{
				out[i] = 0;
				h = Hconn + ncheck*(i - K1 - z);
				l = 0;
				while (h[l]< i && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
				if (h[l] != i)
				{
					printf("Error\n");
				}
			}

			// HARQ section
			for (i = K1 + ccore*z; i< N1; i++)// Core
			{
				out[i] = 0;
				h = Hconn + ncheck*(i - K1);
				l = 0;
				while (h[l]< i && h[l] != -1)
				{
					out[i] ^= out[h[l]];
					l++;
				}
				if (h[l] != i)
				{
					printf("Error\n");
				}
			}
			break;



		}
		if (true)
		{
			bool OK;
			if (period != 0 || shorten > 0)OK=IsCodeword(buff);
			else                           OK=IsCodeword(outt);
//			if (OK)printf("OK\n");
//			else   printf("NOT OK\n");
		}
		if(period!=0 || shorten > 0)
		{
			if (lastsh == 0)
			{
				for(i=shorten,l=0;i<N1;i++)
				{
					if(pattern[i%period]==1) outt[l++]=buff[i];
				}
			}
			else
			{
				l = 0;
				for (i = 0; i < K; i++)
				{
					if(pattern[i%period]==1) outt[l++]=buff[i];
				}
				for (i = K; i < N1-shorten; i++)
				{
					if(pattern[i%period]==1) outt[l++]=buff[i+shorten];
				}
			}
		}
		outt+=N;
		inp +=K;
	}
}
void LDPC_Encoder::Repeat(const int repeat,	int* out)
{
	int i;
	for(i=0;i<repeat;i++)
	{
		out[i]=bufHARQ[order[counter++]];
		if(counter==N)counter=0;
	}
}
void LDPC_Encoder::LLRCombine(const int repeat,int* llr,const int* newllr)
{
	int i;
	counter=(counter-repeat+N)%N; // go back by repeat
	for(i=0;i<repeat;i++)
	{
		llr[order[counter++]]+=newllr[i];
		if(counter==N)counter=0;
	}

}
void LDPC_Encoder::SetHARQ()
{
	HARQ=true;
	// Partition nodes according to highest degrees
	int i,j;
	int* degv=(int*)calloc(N,4);
	for(i=0;i<N-K;i++)
	{
		for(j=0;j<this->ncheck;j++)
		{
			if(Hconn[i*ncheck+j]>=0)degv[Hconn[i*ncheck+j]]++;
			else break;
		}
	}
	delete[] order; order=new int[N];
	delete[] bufHARQ; bufHARQ=new int[N];
	for(i=0;i<N;i++)order[i]=i;
	int temp;
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(degv[order[j]]>degv[order[i]])
			{
				temp=order[i];
				order[i]=order[j];
				order[j]=temp;
			}
		}
	}
	delete[] degv;
}

void LDPC_Encoder::SetVarPerm(const int type)
{
	delete[] permvar; permvar = 0;
	switch (type)
	{
	case 0: //Random
		permvar = Random_Permutation(N1,123598);
		break;
	case 1: //Palindrome
		if (N1 != 0)permvar = new int[N1];
		for (int i = 0; i < N1 / 2; i++)
		{
			permvar[2 * i] = i;
			permvar[2 * i + 1] = N1 - i;
		}
		break;
	case 2: //Row-column
		if (N1 != 0)permvar = new int[N1];
		for (int i = 0; i < N1 / 2; i++)
		{
			permvar[2 * i] = i;
			permvar[2 * i + 1] = i + N1 / 2;
		}
		break;

	case 3: //Identity
		if (N1 != 0)permvar = new int[N1];
		for (int i = 0; i < N1; i++)permvar[i] = i;
		break;

	}
}

/***************************************************************/
/******************** Interface functions **********************/
/***************************************************************/
LDPC_Encoder* LDPC_From_File(char * namefile)
{
	int i,j,k;
	int r;
	int temp;
	LDPC_Encoder* ldpc=new LDPC_Encoder;
	FILE* a;
	fopen_s(&a, namefile, "r");
	if(a==0)
	{
		printf("LDPC_From_File: File %s not found.\n",namefile);
		exit(1);
	}
	char line[132];
	fgets(line, sizeof(line), a);
	fscanf(a,"#\t%d",&ldpc->N1);
	fscanf(a,"\t%d",&r);
	fscanf(a, "\t%d", &ldpc->N);
	fscanf(a, "\t%d", &ldpc->K);
	fscanf(a, "\t%d", &ldpc->z);
	fscanf(a, "\t%d", &ldpc->gap);
	fscanf(a, "\t%d", &ldpc->ncheck);
	fscanf(a,"\n");
	ldpc->K1 = ldpc->N1 - r;
	ldpc->shorten = ldpc->K1-ldpc->K;
	ldpc->Hconn  = new int[ldpc->N1*ldpc->ncheck];

	int*h =ldpc->Hconn;
	int z = ldpc->z;
	int s, c;

	for(i=0;i<r;i+=z)
	{
		for(j=0;j<ldpc->ncheck;j++)
		{
			fscanf_s(a,"%d",&temp);
			if(temp<0)
			{
				for (k = 0; k < z; k++)h[j+k*ldpc->ncheck]=-1;
				break;
			}
			else
			{
				c = (temp / z); 
				s = temp%z;
				for (k = 0; k < ldpc->z; k++)
				{
					h[j + k*ldpc->ncheck] = z*c + (s + k) % z;
				}
			}
		}
		h+=z*ldpc->ncheck;
		fscanf(a,"\n");
	}
	if (ldpc->N1>ldpc->N)
	{
		int period = ldpc->N1 - ldpc->shorten;
		int* pattern = new int[period];
		for (i = 0; i < period; i++)
		{
			fscanf(a, "%1d", &pattern[i]);
			if (i % 100 == 99)fscanf(a, "\n");
		}
		fscanf(a, "\n");
		ldpc->SetPuncturing(period, pattern);
		delete[] pattern;
	}
	if (ldpc->gap > 0)
	{
		unsigned int t2=0;
		delete[]ldpc->gg;
		ldpc->gg = new int[ldpc->K1*ldpc->gap];
		for (i = 0; i < ldpc->gap; i++)
		{
			for (k = 0; k < ldpc->K1; k++)
			{
				if(k%32==0)if(fscanf(a, "%x,",&t2)==0)printf("Error.\n");
				ldpc->gg[i*ldpc->K1 + k]=(t2>>(k%32))&1;
//				fscanf(a, "%1d",&ldpc->gg[i*ldpc->K1 + k]);
			}
			fscanf(a, "\n");
		}
		ldpc->SetEncodingMethod(5);
	}
	fclose(a);
	return ldpc;
}
void BuildGeneratingMatrix(int* h,const int N, const int NK, const int ncheck, int* G, int* p)
{
	int i,l,j,temp,a1,a2;
	int K=N-NK;
	int *H;
	H = (int*)calloc(NK*N,sizeof(int));
	for(i=0;i<N;i++)p[i]=i;

	// Fill Parity check matrix
	for(i=0;i<NK;i++)
	{
		for(l=0;l<ncheck;l++)
		{
			if(h[i*ncheck+l]>=0)H[i*N+h[i*ncheck+l]]=1;
		}
	}

	// Forward reduction
	for(i=0;i<NK;i++)
	{
		// Look for the pivot 
		for(l=i;l<N;l++)
		{
			if(H[i*N+p[l]]!=0)
			{
				temp=p[l];
				p[l]=p[i];
				p[i]=temp;
				break;
			}
		}
		if(l==N)
		{
			continue;
			//printf("%d\t%d\n",i,j);
			//exit(1);
		}

			
		// Zeroes all following rows
		a2 = H[i*N+p[i]];
		printf("%d,%d\n",i,a2);
		for(j=i+1;j<NK;j++)
		{
			if(H[j*N+p[i]]!=0)
			{
//				printf("%d\t%d\n",i,j);
				a1 = H[j*N+p[i]];
				for(l=0;l<N;l++)
				{
					H[j*N+l] = (a2*H[j*N+l]-a1*H[i*N+l]);
				}
			}
		}
	}


	// Backward reduction
	for(i=NK-1;i>0;i--)
	{		
		// Zeroes all previous rows
		a2 = H[i*N+p[i]];
		for(j=i-1;j>=0;j--)
		{
			if(H[j*N+p[i]]!=0)
			{
				a1 = H[j*N+p[i]];
				for(l=0;l<N;l++)
				{
					H[j*N+l] = (a2*H[j*N+l]-a1*H[i*N+l]);
				}
			}
		}
	}
	// Fill G
	for(i=0;i<NK;i++)
	{
		for(j=0;j<K;j++)
		{
			G[i*K+j]=H[i*N+p[j+NK]];
		}
	}
	delete[] H;
}
LDPC_Encoder* LDPC_From_Generating_Matrix(char * namegenmat,char * nameparmat)
{
	int i,j,k;
	int l;
	int K,N,NK,ncheck;
	int *G;


	/* Loading of sparse H */
	FILE* b; fopen_s(&b, nameparmat, "r");
	if(b==0)
	{
		printf("LDPC_From_Generating_Matrix: File %s not found.\n",nameparmat);
		exit(1);
	}
	fscanf_s(b,"%d %d\n",&NK,&N);
	fscanf_s(b,"%d",&ncheck);
	int *h= new int[N*ncheck];
	K=N-NK;

	int temp;
	for(i=l=0;i<NK;i++)
	{
		if(i>0)
		{
			fscanf_s(b,"%d",&temp);
			if(temp!=ncheck)exit(1);
		}
		for(j=0;j<ncheck;j++)
		{
			if(fscanf_s(b,"%d",&h[l])!=1)
			{
				printf("LDPC_From_Generating_Matrix: Reading error in file %s.\n",nameparmat);
				exit(1);
			};

			// Inserted for Ericsson code... (Cyclic shift of columns)
			//if(h[l]<NK)h[l] +=K;
			//else	   h[l]-=NK;
			l++;

		}
		fscanf_s(b,"\n");
	}
	fclose(b);


	/* Loading of Dense systematic G */
	FILE* a; fopen_s(&a, namegenmat, "r");
	if(a==0)
	{
		int* p;
		printf("LDPC_From_Generating_Matrix: File %s not found.\n",namegenmat);
		printf("LDPC_From_Generating_Matrix: Build Generating Matrix...");
		G = (int*)calloc(NK*K,4);
		p = new int[N];
		BuildGeneratingMatrix(h,N,NK,ncheck,G,p);
		printf("done.\n");

		// Sorting columns (first systematic then parity check)
		for(i=0;i<NK*ncheck;i++)
		{
			if(h[i]>=0)
			{
				if(p[h[i]]>=NK)
				{
					h[i]=p[h[i]]-NK;
				}
				else
				{
					h[i]=p[h[i]]+K;
				}
			}
		}
		delete[] p;
	}
	else
	{
		G=new int[K*NK];
		fscanf_s(a,"%d %d\n",&NK,&K);
		for(i=l=0;i<NK;i++)
		{
			for(k=0;k<K;k++)
			{
				if(fscanf_s(a,"%1d",&G[l++])!=1)
				{
					printf("LDPC_From_Generating_Matrix: Reading error in file %s.\n",namegenmat);
					exit(1);
				};
			}
			fscanf_s(a,"\n");
		}
		fclose(a);
	}

	LDPC_Encoder* ldpc=new LDPC_Encoder;
	ldpc->SetEncodingMethod(1);
	ldpc->K     = K;
	ldpc->N	    = ldpc->N1    =N;
	ldpc->Hconn = h;
	ldpc->G     = G;
	ldpc->ncheck = ncheck;
	return ldpc;
}
#include "cppincludes/DVBcodes.cpp"
LDPC_Encoder* LDPC_DVBS2(int type, int shortblock)
{
	int q;
	const int *hh=0;
	LDPC_Encoder* enc=new LDPC_Encoder;
	bool error =false;

	switch(shortblock)
	{
	case 0:enc->N=enc->N1=64800;break;
	case 1:enc->N=enc->N1=16200;break;
	case 2:enc->N=enc->N1=4096;break;
	default: error=true;
	}

	int dmax;
	switch(type)
	{
	case 14:
		printf("LDPC DVBS2: Rate 1/4 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong14;q=135;enc->K=16200;dmax=2;break;
		case 1:hh=HconnDVBshort14 ;q=36 ;enc->K=3240; dmax=2;break;
		case 2:hh=HconnDVBveryshort14;q=24 ;enc->K=1024; dmax=2;break;
		default: error=true;
		}
		break;
	case 13:
		printf("LDPC DVBS2: Rate 1/3 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong13;q=120;enc->K=21600;dmax=3;break;
		case 1:hh=HconnDVBshort13;q=30; enc->K=5400 ;dmax=3;break;
		default: error=true;
		}
		break;
	case 25:
		printf("LDPC DVBS2: Rate 2/5 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong25;q=108;enc->K=25920;dmax=4;break;
		case 1:hh=HconnDVBshort25;q=27; enc->K=6480 ;dmax=4;break;
		default: error=true;
		}
		break;
	case 12:
		printf("LDPC DVBS2: Rate 1/2 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong12;q=90;enc->K=32400;dmax=5;break;
		case 1:hh=HconnDVBshort12;q=25;enc->K=7200 ;dmax=5;break;
		case 2:hh=HconnDVBveryshort12;q=16;enc->K=2048 ;dmax=5;break;
		default: error=true;
		}
		break;
	case 35:
		printf("LDPC DVBS2: Rate 3/5 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong35;q=72;enc->K=38880;dmax=9;break;
		case 1:hh=HconnDVBshort35;q=18;enc->K=9720 ;dmax=9;break;
		default: error=true;
		}
		break;
	case 23:
		printf("LDPC DVBS2: Rate 2/3 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong23;q=60;enc->K=43200;dmax=8;break;
		case 1:hh=HconnDVBshort23;q=15;enc->K=10800;dmax=8;break;
		default: error=true;
		}
		break;
	case 34:
		printf("LDPC DVBS2: Rate 3/4 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong34;q=45;enc->K=48600;dmax=12;break;
		case 1:hh=HconnDVBshort34;q=12;enc->K=11880;dmax=11;break;
		case 2:hh=HconnDVBveryshort34;q=8;enc->K=3072;dmax=15;break;
		default: error=true;
		}
		break;
	case 45:
		printf("LDPC DVBS2: Rate 4/5 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong45;q=36;enc->K=51840;dmax=16;break;
		case 1:hh=HconnDVBshort45;q=10;enc->K=12600;dmax=11;break;
		default: error=true;
		}
		break;
	case 56:
		printf("LDPC DVBS2: Rate 5/6 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong56;q=30;enc->K=54000;dmax=20;break;
		case 1:hh=HconnDVBshort56;q=8 ;enc->K=13320;dmax=17;break;
		default: error=true;
		}
		break;
	case 89:
		printf("LDPC DVBS2: Rate 8/9 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong89;q=20;enc->K=57600;dmax=25;break;
		case 1:hh=HconnDVBshort89;q=5 ;enc->K=14400;dmax=25;break;
		default: error=true;
		}
		break;
	case 910:
		printf("LDPC DVBS2: Rate 9/10 (N=%d)\n",enc->N);
		switch(shortblock)
		{
		case 0:hh=HconnDVBlong910;q= 18;enc->K=58320;dmax=28;break;
		default: error=true;
		}
		break;
	default:
		error=true;
	}

	if(error)
	{
		printf("LDPC_DVBS2: Code not available.\n");
		exit(1);
	}
	enc->K1 = enc->K;
	enc->N1 = enc->N;
	dmax+=2;
	int r = enc->N - enc->K;
	int *cdeg=new int[r];


	enc->ncheck=dmax;
	enc->Hconn=new int[dmax*r];

	int zexp;
	if(shortblock==2)zexp=128;
	else			 zexp=360;
	int p;
	int i,j,par;
	i=0;
	for(i=0;i<r;i++)cdeg[i]=0;
	for(i=0;i<enc->K/zexp;i++)
	{
		while((p=*hh++) != -1)
		{
			if(p<0 || p>=r)
			{
				printf("!!");
			}
			p%=r;
//			fscanf_s(f,"%d%c\n",&p,&c);
			enc->Hconn[p*dmax+cdeg[p]++]=i*zexp;
			for(j=1;j<zexp;j++)
			{
				par=(p+j*q)%r;
				enc->Hconn[par*dmax+cdeg[par]++]=i*zexp+j;
			}
		}
		//while(c!='\n');
	}
	enc->Hconn[cdeg[0]++]=enc->K;
	if(cdeg[0]!=dmax)
		enc->Hconn[cdeg[0]]=-1;
	for(i=1;i<r;i++)
	{
			enc->Hconn[i*dmax+cdeg[i]++] =enc->K+i-1;
			enc->Hconn[i*dmax+cdeg[i]++] =enc->K+i;
			if(cdeg[i]<dmax)
			{
				enc->Hconn[i*dmax+cdeg[i]++] =-1;
			}
			else if(cdeg[i]>dmax)
			{
				printf("Error in LDPC_DVBS2 (%d).\n",cdeg[i]);
				exit(1);
			}
	}
	delete[] cdeg;
	return enc;

}
#include "cppincludes/DVBcodes_SX.cpp"
LDPC_Encoder* LDPC_DVBSX(const int type,const int shortblock)
{
	int q;
	const int *hh=0;
	LDPC_Encoder* enc=new LDPC_Encoder;
	bool error =false;
	int dmax=28;
	int Xs,P,Xp;
	Xs=Xp=P=0;
	switch(shortblock)
	{
	case 0:// Long
		enc->N=enc->N1=64800;
		switch(type)
		{
		case 15:hh=HconnDVBSX_long15;              enc->K=12960;q=144;break;
		case 29:hh=HconnDVBSX_long29;              enc->K=14400;q=140;Xs=0;P=15;Xp=3240;break;
		case 1345:hh=HconnDVBSX_long1345;          enc->K=18720;q=128;break;
		case 920:hh=HconnDVBSX_long920;            enc->K=29160;q=99; break;
		case 90180:hh=HconnDVBSX_long90180;        enc->K=32400;q=90; break;
		case 96180:hh=HconnDVBSX_long96180;        enc->K=34560;q=84; break;
		case 1120:hh=HconnDVBSX_long1120;          enc->K=35640;q=81; break;
		case 100180:hh=HconnDVBSX_long100180;      enc->K=36000;q=80; break;
		case 104180:hh=HconnDVBSX_long104180;      enc->K=37440;q=76; break;
		case 2645:hh=HconnDVBSX_long2645;          enc->K=37440;q=76; break;
		case 1830:hh=HconnDVBSX_long1830;			enc->K=38880;q=72; break;
		case 2845:hh=HconnDVBSX_long2845;          enc->K=40320;q=68; break;
		case 2336:hh=HconnDVBSX_long2336;          enc->K=41400;q=65; break;
		case 116180:hh=HconnDVBSX_long116180;      enc->K=41760;q=64; break;
		case 2030:hh=HconnDVBSX_long2030;			enc->K=43200;q=60; break;
		case 124180:hh=HconnDVBSX_long124180;      enc->K=44640;q=56; break;
		case 2536:hh=HconnDVBSX_long2536;          enc->K=45000;q=55; break;
		case 128180:hh=HconnDVBSX_long128180;      enc->K=46080;q=52; break;
		case 1318:hh=HconnDVBSX_long1318;          enc->K=46800;q=50; break;
		case 132180:hh=HconnDVBSX_long132180;      enc->K=47520;q=48; break;
		case 2230:hh=HconnDVBSX_long2230;			enc->K=47520;q=48; break;
		case 135180:hh=HconnDVBSX_long135180;      enc->K=48600;q=45; break;
		case 79:hh=HconnDVBSX_long79;              enc->K=50400;q=40; break;
		case 140180:hh=HconnDVBSX_long140180;      enc->K=50400;q=40; break;
		case 154180:hh=HconnDVBSX_long154180;      enc->K=55440;q=26; break;
		case 144180:hh=HconnDVBSX_long144180;      enc->K=51840;q=36; break;
		case 150180:hh=HconnDVBSX_long150180;      enc->K=54000;q=30; break;
		default: error=true;
		}
		break;
	case 1:// Short
		enc->N=enc->N1=16200;
		switch(type)
		{
		case 15:hh=HconnDVBshort14 ;q=36 ;enc->K=3240; dmax=2;Xs=560;P=30;Xp=250;break;
		case 1145 :hh=HconnDVBSX_short1145;enc->K=3960 ;q=34;break;
		case 11451:hh=HconnDVBSX_short1145;enc->K=3960 ;q=34;;Xs=0;P=15;Xp=810;break; //For 11/45  
		case 415:hh=HconnDVBSX_short415;enc->K=4320   ;q=33;break;
		case 1445:hh=HconnDVBSX_short1445;enc->K=5040 ;q=31;break;
		case 715:hh=HconnDVBSX_short715;enc->K=7560   ;q=24;break;
		case 815:hh=HconnDVBSX_short815;enc->K=8640   ;q=21;break;
		case 2645:hh=HconnDVBSX_short2645;enc->K=9360 ;q=19;break;
		case 3245:hh=HconnDVBSX_short3245;enc->K=11520;q=13;break;
		default: error=true;
		}
		break;
	case 2:// medium
		enc->N=enc->N1=32400;// Puncturing
		switch(type)
		{
		case 15 :hh=HconnDVBSX_medium15;enc->K=5840+640 ;q=72;	Xs=640;P=25;Xp=980;break;
		case 1145:hh=HconnDVBSX_medium1145;enc->K=7920   ;q=68;	Xs=0  ;P=15;Xp=1620;break;
		case 13:hh=HconnDVBSX_medium13;enc->K=10800 ;q=60;		Xs=0;  P=13;Xp=1620;break;
		default: error=true;
		}
		break;
	default: error=true;
	}
	enc->K1 = enc->K;

	if(error)
	{
		printf("LDPC_DVBSX Code not available.\n");
		delete[] enc;
		enc=0;
		return enc;
	}

	dmax+=2;
	int r = enc->N - enc->K;
	int *cdeg=new int[r];


	enc->ncheck=dmax;
	enc->Hconn=new int[dmax*r];

	int zexp;
	zexp=360;
	if(enc->K %zexp != 0)
	{
		printf("Error in module LDPC_DVBSX\n");
		exit(1);
	}
	int p;
	int i,j,par;
	i=0;
	for(i=0;i<r;i++)cdeg[i]=0;
	for(i=0;i<enc->K/zexp;i++)
	{
		while((p=*hh++) != -1)
		{
			//if(p<=0|| p>=r)
			//{
			//	printf("!!");
			//}
			p %= r;

			enc->Hconn[p*dmax + cdeg[p]++]=i*zexp;
			for(j=1;j<zexp;j++)
			{
				par=(p+j*q)%r;
				enc->Hconn[par*dmax+cdeg[par]++]=i*zexp+j;
			}
		}
	}
	enc->Hconn[cdeg[0]++]=enc->K;
	if(cdeg[0]!=dmax)
		enc->Hconn[cdeg[0]]=-1;
	for(i=1;i<r;i++)
	{
			enc->Hconn[i*dmax+cdeg[i]++] =enc->K+i-1;
			enc->Hconn[i*dmax+cdeg[i]++] =enc->K+i;
			if(cdeg[i]<dmax)
				enc->Hconn[i*dmax+cdeg[i]++] =-1;
			else if(cdeg[i]>dmax)
			{
				printf("Error in LDPC_DVBSX (%d).\n",cdeg[i]);
				exit(1);
			}
	}
	int* pattern;
	if(Xp>0)
	{
		pattern=new int[enc->N1];
		int i,l;
		for(i=0;i<enc->N1;i++)pattern[i]=1;

		for(i=enc->K,l=0;l<Xp;l++,i+=P)
		{
			pattern[i]=0;
		}
		enc->SetPuncturing(enc->N1,pattern);
		delete[] pattern;
		if(Xs>0)// Shortening
		{
			enc->Shorten(Xs);
		}
	}
	delete[] cdeg;
	return enc;
}
LDPC_Encoder* LDPC_RegularRandom2(const int K,const int N,int ncheck)
{
	int i,j,ipick,seed,temp;

	seed=712991519;

	if(ncheck<=0)ncheck=(int)((double)(3.*N)/(N-K)+0.9999);
	int* Hconn=new int[ncheck*(N-K)];
	int *vect=new int[N];
	int l,save;
	l=0;
	for(i=0;i<N;i++)vect[i]=i;
	for(i=0;i<N-K;i++)
	{
		for(j=0;j<ncheck;j++)
		{
			if(j==0){ipick =(unif_int(seed)%(N-l-ncheck));save=ipick;}
			else	
			{
				do
				{
					ipick = (save%ncheck) + ncheck*(unif_int(seed)%((N-l-ncheck)/ncheck));
				}while(ipick==save);
			}
			// remove short loops
			ipick+=l+ncheck;
			temp=vect[l+j];
			vect[l+j]=vect[ipick];
			vect[ipick]=temp;
			Hconn[i*ncheck+j]=vect[l+j];
//			l=(l+1)%N;
		}
		l=(l+ncheck)%N;
		/* Sort [optional]
		for(j=0;j<ncheck-1;j++)
		{
			for(jj=j+1;jj<ncheck-1;jj++)
			{
				if(Hconn[i*ncheck+j]>Hconn[i*ncheck+jj])
				{
					temp=Hconn[i*ncheck+jj];
					Hconn[i*ncheck+jj]=Hconn[i*ncheck+j];
					Hconn[i*ncheck+j]=temp;
				}
			}
		} */
	}
	LDPC_Encoder *a=new LDPC_Encoder;
	a->SetParameters(K,N,ncheck,Hconn);
	delete[] Hconn;
	delete[] vect;
	return a;
}
LDPC_Encoder * LDPC_NRQUALCOMM(const int K, const int N, int &family, const int sol)
{
	#include "cppincludes/NRLTECodes.cpp"
	int kmin, kmax, cmin, cmax,ccor;
	const int *fam;
	const int *conn=0;
	double rate = (double)K / N;
	if (family <= 0)
	{
		if  (rate<= 2. / 5.)family = 1;     // 0.4
		else if (rate<= 2. / 3.)family = 2; // 0.666
		else if (rate<= 1.)family = 3;		// 8/9
	}
	if (K > 8960 && family == 1)
	{
		printf("Warning: forcing to wrong family 2 for size problems.\n");
		family = 2;// force to second family
	}
	if (K > 17920 && family == 2)
	{
		printf("Warning: forcing to wrong family 3 for size problems.\n");
		family = 3;// force to third family
	}
	if (K > 26880)
	{
		printf("Warning: Too large block size.\n");
		family = 0;
	}

	if (family == 0)return 0;
	fam = 0;
	switch (family)
	{
	case 1:fam = Low;     conn = LowP;break;
	case 2:fam = Middle;  conn = MiddleP;break;
	case 3:fam = High;    conn = HighP; break;
	}
	kmax = fam[0];
	kmin = fam[1];
	ccor = fam[2];
	cmax = fam[3];
	cmin = fam[4];

	//Find all solutions Z_i∈{ 8,…,896 } so that k_(b, min)≤K / Z_i ≤k_(b, max)
	int Z;
	int Zmax = K / kmin; // Zmin
	int Zmin =(int)((double)K / kmax+0.99999); //Zmax
	int jcl = 1;
	int i;
	int kb;
	int cb, bsol;
	unsigned int S,Sh, P,Zopt,best;
	
	// Find best solution or required solution
	best = 10000;
	int nsol = 0;
	for (jcl = 1; jcl <= 7; jcl++)
	{
		for (i = 4; i <= 7; i++)
		{
			Z = (1 << jcl)*i;
			if (Z > Zmax)goto end;
			if (Z >= Zmin)
			{
//				printf("Z=%d\n", Z);
				//Set the number of base graph information variable nodes to  k_b = ⌈K / Z_i ⌉ by deleting the last k_(b, max) - k_b base information variable nodes from the base graph of the family
				kb = (int)(((double)K / Z) + 0.9999);
				//Append the first ⌈(N - K) / Z_i ⌉ + p_b parity variable nodes unless ⌈(N - K) / Z_i ⌉ + p_b<c_(b, core) in which case the c_(b, core) parity variable nodes are appended. (The number of base check nodes is equal to the number of base parity variable nodes.)
				cb = 2 + (int)((double)(N - K) / Z + 0.9999999);
				if (cb < ccor)cb = ccor;
				S = kb*Z - K;
				P = cb*Z - (N - K + 2 * Z);
//				printf("K=%d\tN=%d\tZ=%d\t kb=%d\t cb=%d, K1=%d\t N1=%d\tS=%d\tP=%d\n", K, N, Z, kb, cb, kb*Z, (kb + cb)*Z, S, P);
				if ((sol==nsol) || (sol==-1) && (S + P < best))
				{
					bsol = nsol;
					best = S + P;
					Zopt = Z;
				}
				nsol++;
			}
		}

	}
end:
	if (best == 10000)return 0;

	Z = Zopt;
	jcl = 0;
	while ((Z >> jcl)>7)jcl++;
	kb = (int)(((double)K / Z) + 0.9999);
	cb = 2 + (int)((double)(N - K) / Z + 0.9999999);
	if (cb < ccor)cb = ccor;
	if (cb > cmax)return 0;
	S = kb*Z - K;
	P = cb*Z - (N - K + 2 * Z);
	//printf("K=%d\tN=%d\nfamily=%d\tZ=%d=%dx%d\nkb=%d\t cb=%d\nK1=%d\t N1=%d\nS=%d\tP=%d (sol=%d/%d)\n", 
	//	K, N, family,Z,1<<jcl,Z/(1<<jcl), kb, cb, kb*Z, (kb + cb)*Z, S, P,bsol,nsol);

	// Upload parameters
	int deg,l;

	int maxdeg = 25;
	int *hconn = new int[maxdeg*cb];
	int *shift = new int[maxdeg*cb];
	l = 1;
	for (i = 0; i < cb; i++)
	{
		deg = 0;
//		printf("l=%d\t", conn[l - 1]);
		while(true)
		{
			if (conn[l] == -1)
			{
				hconn[i*maxdeg + deg] = -1;
				shift[i*maxdeg + deg] = -1;
//				printf("deg=%d\n",deg);
				break;
			}
			else
			{
				if (conn[l] <= kb || conn[l]> kmax)
				{
					if(conn[l]<=kb)hconn[i*maxdeg + deg] = conn[l]-1;
					else           hconn[i*maxdeg + deg] = conn[l]-(kmax-kb)-1;
					Sh=conn[l+1];
	//				Sh = 0x19EF51;jcl = 3;  Example
					shift[i*maxdeg + deg] = ((Sh>>(21-jcl))<<2) + ((Sh>>(2*(jcl-1)))&3);
					deg++;
				}
				l += 2;
			}
		}
		l+=2;
	}

	LDPC_Encoder* ldpc = new LDPC_Encoder;
	if (family == 1)
	{
		if (jcl > 1)// cluster
		{
			ldpc->hhh = 1 << (jcl - 2);
			switch (Z >> jcl)
			{
			case 4: ldpc->aaa = 4  * ldpc->hhh; break;
			case 5: ldpc->aaa = 8  * ldpc->hhh; break;
			case 6: ldpc->aaa = 20 * ldpc->hhh; break;
			case 7: ldpc->aaa = 0  * ldpc->hhh; break;
			}
		}
		else
		{
			ldpc->hhh = 0;
			ldpc->aaa = 4;
		}
		fprintf(stdout, "a=%d,h=%d\n",ldpc->aaa, ldpc->hhh);

	}
	else
	{
		ldpc->aaa = 0;
		ldpc->hhh = 0;
	}
	int* h = new int[cb*Z*maxdeg];
	int* hh;
	int j;
	l = 0;
	for (i = 0; i<cb; i++)
	{
		for (j = 0; j<Z; j++)
		{
			hh = h + maxdeg*(i*Z + j);
			deg = 0;
			while (hconn[i*maxdeg + deg]>=0)
			{
				hh[deg++] = hconn[i*maxdeg + deg]*Z + (shift[i*maxdeg + deg] + j) % Z;
			}
			hh[deg] = -1;
		}
	}
	ldpc->x = 0;  // Exception 
	ldpc->z = Z;
	ldpc->ccore = ccor;
	ldpc->K = ldpc->K1 = kb*Z;
	ldpc->N = ldpc->N1 = (kb+cb)*Z;
	ldpc->Hconn = h;
	ldpc->ncheck = maxdeg;
	ldpc->SetEncodingMethod(4);

	ldpc->Shorten(S,1); // Shortening at the end of data block

	int*pattern = new int[ldpc->N];
	for (i = 0; i < ldpc->N; i++)pattern[i] = 1;
	for (i = 0; i < 2 * Z; i++)pattern[i] = 0;
	for (i = ldpc->N - P; i < ldpc->N; i++)pattern[i] = 0;
	ldpc->SetPuncturing(ldpc->N, pattern);
	delete[] hconn;
	delete[] shift;
	delete[] pattern;
	return ldpc;
}
LDPC_Encoder * LDPC_NR5G(const int K, const int N, int &family, int &cluster)
{
#include "cppincludes/NRLTECodes.cpp"
	int  cbmax,cbcor,kb,kmax;
	const int *conn;
	double rate = (double)K / N;
	bool outofspec = false;
	if (K <= 0)return 0;
	if (N <= 0)return 0;
	if (family != 1 &&family !=2 )
	{
		family = 1;
		if (K<=292)family = 2;                        // BG2
		if ((rate <= 0.67) && (K <= 3824))family = 2;    // BG2
		if ((rate <= 0.67) && (K <= 3840))family = 2;    // BG2
		if ((rate < 0.25))family = 2;    // BG2
	}

	switch (family)
	{
	case 1:
		if (K > 8448)
		{
			printf("Warning code out of specs (K=%d).\n",K);
			outofspec = true;
		}
		kmax = 22;  cbcor = 4;  cbmax = 46; 
		kb = 22;
		break;
	case 2:kmax = 10;  cbcor = 4; cbmax = 42;
		if (K > 3840)
		{
			printf("Warning code out of specs (K=%d).\n",K);
			outofspec = true;
		}
		if (K > 640)	 kb = 10;
		else if (K > 560)kb = 9;
		else if (K > 192)kb = 8;
		else             kb = 6;
		break;
	}

	/*
	find the minimum value of    in all sets of lifting sizes in Table 5.3.2-1, denoted as  , 
	such that  , and denote   for LDPC base graph 1 and   for LDPC base graph 2;
	*/
retry:
	//Find the solution for Z with largest all solutions Z_i∈{ 8,…,896 } so that k_(b, min)≤K / Z_i ≤k_(b, max)
	int Z;
	int jcl = 1;
	int i;
	int cb,iopt;
	unsigned int S, P;

	// Find the parallelism
	int Zc = 100000;
	for (i = 0; i < 8; i++)
	{
		if (cluster >= 0)i = cluster;
		for (jcl = 0; jcl < 15; jcl++)
		{
			Z = shifts5G[i]<<jcl;
//			if (Z > 384)continue;
			if (kb*Z < K)continue;
			//kb = (int)((double)K / Z+0.9999);
			//if (kb > kmax)continue;
			if (Z < Zc)
			{
				iopt = i;
				Zc = Z;
//				kb = (int)((double)K / Z+0.9999);
			}
			break;
		}
		if (cluster >=0)break;
	}
	if (Zc > 384)
	{
		printf("Warning code out of specs (Z=%d).\n",Zc);
		outofspec = true;
	}
	if (Zc == 100000)return 0;
	Z = Zc;
	cb = 2 + (int)((double)(N - K) / Z + 0.9999999);  // Number of columns
	if (cb < cbcor)return 0;
	//{
	//	printf("Warning\n");
	//	cb = cbcor;
	//}
	if (cb > cbmax)
	{
		printf("Warning code out of specs (cb=%d).\n", cb);
		//Zc *= 2;
		//kb = (int)((double)K / Zc+0.9999);
		kb--;
		printf("Decreasing kb to %d.\n", kb);
		outofspec = true;
		goto retry;
	}
	S = kb*Z - K;
	P = cb*Z - (N - K + 2 * Z);
	printf("K=%d\tN=%d\nfamily=%d\tZ=%d=%dx%d\nkb=%d/%d\t cb=%d/%d\nK1=%d\t N1=%d\nS=%d\tP=%d\n",
		K, N, family, Z,  Z / shifts5G[iopt], shifts5G[iopt],kb, kmax,cb,cbmax, kb*Z, (kb + cb)*Z, S, P);
	cluster = iopt;
	//FILE* xxx = fopen("codes.txt", "a");
	//fprintf(xxx,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",
	//	K, N, family, Z,  Z / shifts[iopt], shifts[iopt],kb, cb, kb*Z, (kb + cb)*Z, S, P);
	//fclose(xxx);

	// Pick the correct matrix
	conn = 0;
	switch (family)
	{
	case 1:
		switch (iopt)
		{
		case 0: conn = BG11; break;
		case 1:conn = BG12; break;
		case 2:conn = BG13; break;
		case 3:conn = BG14; break;
		case 4:conn = BG15; break;
		case 5:conn = BG16; break;
		case 6:conn = BG17; break;
		case 7:conn = BG18; break;
		}break;
	case 2:
		switch (iopt)
		{
		case 0:conn = BG21; break;
		case 1:conn = BG22; break;
		case 2:conn = BG23; break;
		case 3:conn = BG24; break;
		case 4:conn = BG25; break;
		case 5:conn = BG26; break;
		case 6:conn = BG27; break;
		case 7:conn = BG28; break;
		}break;
	}
//	int* Hnew=Stronger(22, 26, conn, 1, 4);

	// Construct parity check matrix
	int deg, l,j;

	int maxdeg = kmax+cbcor;
	int *hconn = new int[maxdeg*cb];
	int *shift = new int[maxdeg*cb];
	l = 1;
	for (i = 0; i < cb; i++)
	{
		deg = 0;
		for (j = 0; j < kb; j++)
		{
			if (conn[j] == -1)continue;
			hconn[i*maxdeg + deg] = j;
			shift[i*maxdeg + deg] = (conn[j]%Zc);
			deg++;
		}
		for (j = kmax; j < maxdeg; j++)
		{
			if (conn[j] == -1)continue;
			hconn[i*maxdeg + deg] = j-(kmax-kb);
			shift[i*maxdeg + deg] = (conn[j]%Zc);
			deg++;

		}
		if (i >= cbcor)
		{
			hconn[i*maxdeg + deg] = kb+i;
			shift[i*maxdeg + deg] = 0;
			deg++;
		}
		if (deg < maxdeg) {
			hconn[i*maxdeg + deg] = -1; shift[i*maxdeg + deg] = -1;
		}
		conn += maxdeg;

	}

	LDPC_Encoder* ldpc = new LDPC_Encoder;

	// For encoding the first Z parity check bits
	ldpc->aaa = 0;
	switch (family)
	{
	case 1:
		switch (iopt)
		{
		case 6:ldpc->aaa = 105%Z; break;
		}
		break;
	case 2:
		switch (iopt)
		{
		case 0:case 1:case 2:case 4:case 5:case 6:ldpc->aaa = 1; break;
		}
		break;
	}
	int* h = new int[cb*Z*maxdeg];
	int* hh;
	l = 0;
	for (i = 0; i<cb; i++)
	{
		for (j = 0; j<Z; j++)
		{
			hh = h + maxdeg*(i*Z + j);
			deg = 0;
			while (hconn[i*maxdeg + deg] >= 0)
			{
				hh[deg++] = hconn[i*maxdeg + deg] * Z + (shift[i*maxdeg + deg] + j) % Z;
			}
			hh[deg] = -1;
		}
	}
	ldpc->x = 0;  // Exception 
	ldpc->z = Z;
	ldpc->ccore = cbcor;
	ldpc->K = ldpc->K1 = kb*Z;
	ldpc->N = ldpc->N1 = (kb + cb)*Z;
	ldpc->Hconn = h;
	ldpc->ncheck = maxdeg;
	ldpc->SetEncodingMethod(6); // Set the proper encoding method

	ldpc->Shorten(S, 1); // Shortening at the end of data block

	int*pattern = new int[ldpc->N];
	for (i = 0; i < ldpc->N; i++)pattern[i] = 1;
	for (i = 0; i < 2 * Z; i++)pattern[i] = 0;
	for (i = ldpc->N - P; i < ldpc->N; i++)pattern[i] = 0;
	ldpc->SetPuncturing(ldpc->N, pattern);
	delete[] hconn;
	delete[] shift;
	delete[] pattern;
	if (outofspec)family += 2;
	return ldpc;
}
LDPC_Encoder * LDPC_NR5G_strong(const int K, const int N, int &family, int &cluster)
{
#include "cppincludes/NRLTECodes.cpp"
	int  cbmax, cbcor, kb, kmax;
	const int *conn;
	double rate = (double)K / N;
	bool outofspec = false;
	if (K <= 0)return 0;
	if (N <= 0)return 0;
	int st = 0;// Added line in core matrix for strong code
	if (family != 1 && family != 2)
	{
		family = 1;
		if (K <= 292)family = 2;                        // BG2
		if ((rate <= 0.67) && (K <= 3824))family = 2;    // BG2
		if ((rate <= 0.67) && (K <= 3840))family = 2;    // BG2
		if ((rate < 0.25))family = 2;    // BG2
	}

	switch (family)
	{
	case 1:
		if (K > 8448)
		{
			printf("Warning code out of specs (K=%d).\n", K);
			outofspec = true;
		}
		kmax = 22;  cbcor = 4+st;  cbmax = 46+st;
		kb = 22;
		break;
	case 2:kmax = 10;  cbcor = 4+st; cbmax = 42+st;
		if (K > 3840)
		{
			printf("Warning code out of specs (K=%d).\n", K);
			outofspec = true;
		}
		if (K > 640)	 kb = 10;
		else if (K > 560)kb = 9;
		else if (K > 192)kb = 8;
		else             kb = 6;
		break;
	}



	/*
	find the minimum value of    in all sets of lifting sizes in Table 5.3.2-1, denoted as  ,
	such that  , and denote   for LDPC base graph 1 and   for LDPC base graph 2;
	*/
retry:
	//Find the solution for Z with largest all solutions Z_i∈{ 8,…,896 } so that k_(b, min)≤K / Z_i ≤k_(b, max)
	int Z;
	int jcl = 1;
	int i;
	int cb, iopt;
	unsigned int S, P;

	// Find the parallelism
	int Zc = 100000;
	for (i = 0; i < 8; i++)
	{
		if (cluster >= 0)i = cluster;
		for (jcl = 0; jcl < 15; jcl++)
		{
			Z = shifts5G[i] << jcl;
			//			if (Z > 384)continue;
			if (kb*Z < K)continue;
			//kb = (int)((double)K / Z+0.9999);
			//if (kb > kmax)continue;
			if (Z < Zc)
			{
				iopt = i;
				Zc = Z;
				//				kb = (int)((double)K / Z+0.9999);
			}
			break;
		}
		if (cluster >= 0)break;
	}
	if (Zc > 384)
	{
		printf("Warning code out of specs (Z=%d).\n", Zc);
		outofspec = true;
	}
	if (Zc == 100000)return 0;
	Z = Zc;
	cb = 2 + (int)((double)(N - K) / Z + 0.9999999);  // Number of columns
	if (cb < cbcor)return 0;
	//{
	//	printf("Warning\n");
	//	cb = cbcor;
	//}
	if (cb > cbmax)
	{
		printf("Warning code out of specs (cb=%d).\n", cb);
		//Zc *= 2;
		//kb = (int)((double)K / Zc+0.9999);
		kb--;
		printf("Decreasing kb to %d.\n", kb);
		outofspec = true;
		goto retry;
	}
	S = kb*Z - K;
	P = cb*Z - (N - K + 2 * Z);
	printf("K=%d\tN=%d\nfamily=%d\tZ=%d=%dx%d\nkb=%d/%d\t cb=%d/%d\nK1=%d\t N1=%d\nS=%d\tP=%d\n",
		K, N, family, Z, Z / shifts5G[iopt], shifts5G[iopt], kb, kmax, cb, cbmax, kb*Z, (kb + cb)*Z, S, P);
	cluster = iopt;
	//FILE* xxx = fopen("codes.txt", "a");
	//fprintf(xxx,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",
	//	K, N, family, Z,  Z / shifts[iopt], shifts[iopt],kb, cb, kb*Z, (kb + cb)*Z, S, P);
	//fclose(xxx);

	// Pick the correct matrix
	conn = 0;
	switch (family)
	{
	case 1:
		switch (iopt)
		{
		case 0: conn = BG11; break;
		case 1:conn = BG12; break;
		case 2:conn = BG13; break;
		case 3:conn = BG14; break;
		case 4:conn = BG15; break;
		case 5:conn = BG16; break;
		case 6:conn = BG17; break;
		case 7:conn = BG18; break;
		}break;
	case 2:
		switch (iopt)
		{
		case 0:conn = BG21; break;
		case 1:conn = BG22; break;
		case 2:conn = BG23; break;
		case 3:conn = BG24; break;
		case 4:conn = BG25; break;
		case 5:conn = BG26; break;
		case 6:conn = BG27; break;
		case 7:conn = BG28; break;
		}break;
	}


	// Construct parity check matrix
	int deg, l, j;

	int maxdeg = kmax + cbcor;
	int *hconn = new int[maxdeg*cb];
	int *shift = new int[maxdeg*cb];
	l = 1;
	for (i = 0; i < cb; i++)
	{
		deg = 0;
		if (i <cbcor-st-1 || i>=cbcor)
		{
			for (j = 0; j < kb; j++)
			{
				if (conn[j] == -1)continue;
				hconn[i*maxdeg + deg] = j;
				shift[i*maxdeg + deg] = (conn[j] % Zc);
				deg++;
			}
			for (j = kmax; j < kmax + 1+st ; j++)// Strong lines
			{
				if (j == kmax)
				{
					if (conn[j-st] == -1)continue;
					hconn[i*maxdeg + deg] = j - (kmax - kb);
					shift[i*maxdeg + deg] = (conn[j-st] % Zc);
					deg++;
				}
				else
				{
				}
			}
			for (j = kmax+1+st; j < maxdeg; j++)
			{
				if (conn[j-st] == -1)continue;
				hconn[i*maxdeg + deg] = j - (kmax - kb);
				shift[i*maxdeg + deg] = (conn[j-st] % Zc);
				deg++;
			}
			if (i >= cbcor)// HARQ section
			{
				hconn[i*maxdeg + deg] = kb + i;
				shift[i*maxdeg + deg] = 0;
				deg++;
			}
			conn += maxdeg;
		}
		else
		{
			int sel = (i - (cbcor - st - 1)) % (st + 1);
			int dd;
			for (j = dd = 0; j < kb; j++)
			{
				if (conn[j] == -1)continue;
				if (dd++%(st+1) != sel)continue;
				hconn[i*maxdeg + deg] = j;
				shift[i*maxdeg + deg] = (conn[j] % Zc);
				deg++;
			}
			for (j = kmax; j < kmax + 1 + st; j++)// Strong lines
			{
				if (j == kmax)
				{
					if (conn[j - st] == -1)continue;
					hconn[i*maxdeg + deg] = j - (kmax - kb);
					shift[i*maxdeg + deg] = (conn[j - st] % Zc);
					deg++;
				}
				else
				{
				}
			}
			for (j = kmax + 1 + st; j < maxdeg; j++)// Double diagonal
			{
				if (conn[j - st] == -1)continue;
				hconn[i*maxdeg + deg] = j - (kmax - kb);
				shift[i*maxdeg + deg] = (conn[j - st] % Zc);
				deg++;
			}
		}
		if (deg < maxdeg) {
			hconn[i*maxdeg + deg] = -1; shift[i*maxdeg + deg] = -1;
		}

	}

	LDPC_Encoder* ldpc = new LDPC_Encoder;

	// For encoding the first Z parity check bits
	ldpc->aaa = 0;
	switch (family)
	{
	case 1:
		switch (iopt)
		{
		case 6:ldpc->aaa = 105 % Z; break;
		}
		break;
	case 2:
		switch (iopt)
		{
		case 0:case 1:case 2:case 4:case 5:case 6:ldpc->aaa = 1; break;
		}
		break;
	}
	int* h = new int[cb*Z*maxdeg];
	int* hh;
	l = 0;
	for (i = 0; i<cb; i++)
	{
		for (j = 0; j<Z; j++)
		{
			hh = h + maxdeg*(i*Z + j);
			deg = 0;
			while (hconn[i*maxdeg + deg] >= 0)
			{
				hh[deg++] = hconn[i*maxdeg + deg] * Z + (shift[i*maxdeg + deg] + j) % Z;
			}
			hh[deg] = -1;
		}
	}
	ldpc->x = 0;  // Exception 
	ldpc->z = Z;
	ldpc->ccore = cbcor;
	ldpc->K = ldpc->K1 = kb*Z;
	ldpc->N = ldpc->N1 = (kb + cb)*Z;
	ldpc->Hconn = h;
	ldpc->ncheck = maxdeg;
	ldpc->SetEncodingMethod(6); // Set the proper encoding method

	ldpc->Shorten(S, 1); // Shortening at the end of data block

	int*pattern = new int[ldpc->N];
	for (i = 0; i < ldpc->N; i++)pattern[i] = 1;
	for (i = 0; i < 2 * Z; i++)pattern[i] = 0;
	for (i = ldpc->N - P; i < ldpc->N; i++)pattern[i] = 0;
	ldpc->SetPuncturing(ldpc->N, pattern);
	delete[] hconn;
	delete[] shift;
	delete[] pattern;
	if (outofspec)family += 2;
	return ldpc;
}


int* Stronger(const int k, const int n, const int* H, const int st1, const int st2)
{
	//extension matrices
	int X1[] = {
		1,
		1,
		1,
		0 };
	int X2[] = {
		0, 0,
		1, 0,
		1, 1,
		1, 1,
		0, 1 };
	int X3[] = {
		1, 1, 0,
		0, 0, 0,
		0, 1, 1,
		0, 0, 1,
		1, 0, 0,
		1, 1, 1 };
	int X4[] = {
		0, 1, 0, 0,
		1, 0, 1, 1,
		0, 1, 0, 0,
		0, 0, 1, 1,
		0, 0, 0, 0,
		1, 1, 0, 1,
		1, 0, 1, 0};
	int X5[] = {
		0, 1, 0, 0, 1,
		1, 1, 1, 1, 1,
		0, 0, 0, 0, 0,
		0, 0, 1, 1, 0,
		0, 0, 0, 0, 1,
		1, 0, 0, 1, 0,
		1, 0, 1, 0, 0,
		0, 1, 0, 0, 0 };

	//Add st rows and columns to increase the average column degree
	int n1 = n + (st2 - st1);
	int m1 = (n - k + st2 - st1);
	int mm = (n - k - st1);
	int* Hnew = new int[m1*n1];
	int i, j;
	for (i = 0; i < m1*n1; i++)Hnew[i] = -1;
	for (i = 0; i < n - k-st1; i++)// First rows
	{
		for (j = 0; j < k + st1; j++)Hnew[n1*i + j] = H[n*i + j];
		for (; j < n; j++)	 Hnew[n1*i + j+(st2-st1)] = H[n*i + j];
	}
	/* spread ones in systematic part of
	last st1 rows on all st2 rows */
	int dd = 0;
	int row;
	for (; i < n - k; i++)
	{
		for (j = 0; j < k; j++)
		{
			if(H[n*i + j]<0)continue;
			row = (dd++%st2);
			Hnew[(n - k - st1 + row)*n1 + j] = H[n*i + j];
		}
		for (j = k ; j < k + st1; j++)
			Hnew[n1*i + j] = H[n*i + j];
		for (j = k + st1; j < n; j++)
			Hnew[n1*i + j + (st2 - st1)] = H[n*i + j];
	}
	// Now Insert the BD matrices
	for (i=st1;i<st2;i++)
	{
		Hnew[n1*(n-k+i-st1)+ (k+st2)+(n-k-st1-1+i)%(n-k-st1)] = 0;
	}
	


	return Hnew;
}
#include "cppincludes/WiFiCodes.cpp"
LDPC_Encoder* LDPC_WiFi(const int rate,const int N	)
{
	int z;
	int* H,*tab;
	int K;
	switch(N)
	{
		case 648:
			H=WiFi_648EncTable;
			break;

		case 1296:
			H=WiFi_1296EncTable;
			break;

		case 1944:
			H=WiFi_1944EncTable;
			break;

		default:
			printf("Warning:  %d is Not an 802.11ac LDPC code length.\n",N);
			if(N%24!=0)exit(1);
			if(N<=1296)		H=WiFi_648EncTable;
			else if(N<1944)	H=WiFi_1296EncTable;
			else			H=WiFi_1944EncTable;
			break;

	}
	int ncheck=24;
	z=N/24;
	switch(rate)
	{

	case 12:
		K=N/2;
		tab = H;
		break;

	case 23:
		K=(N*2)/3;
		tab =H+ 12*24;
		break;

	case 34:
		K=(N*3)/4;
		tab=H+(12+8)*24;
		break;

	case 56:
		K=(N*5)/6;
		tab=H+(12+8+6)*24;
		break;
	case 78:
		printf("Warning: 78 is not an 802.11ac LDPC code rate.\n");
		K=(N*7)/8;
		tab=H+(12+8+6+4)*24; // Skip the first line
		break;
	default:
		printf("LDPC_WiFi: Rate %d not available.\n",rate);
		exit(1);
	}

	int r=(N-K)/z;
	int* h=new int[(N-K)*ncheck];
	int*hh;
	int i,j,k,l;
	int shift;
	l=0;
	for(i=0;i<r;i++)
	{
		for(j=0;j<z;j++)
		{
			hh=h + ncheck*(i*z+j);
			l=0;
			for(k=0;k<24;k++)
			{
				if(tab[i*24+k]>=0)
				{
					shift=(int)(((double)tab[i*24+k]*z)/96.);
					hh[l++] = k*z+ (shift+j)%z;
				}

			}
			hh[l]=-1;
		}
		
	}
	int x = 0;
	LDPC_Encoder* ldpc=new LDPC_Encoder;
	ldpc->x     = x;
	ldpc->z     = z;
	ldpc->K     = K;
	ldpc->N	    = ldpc->N1    =N;
	ldpc->Hconn = h;
	ldpc->ncheck = ncheck;
	ldpc->SetEncodingMethod(2);
	return ldpc;
}
int liftx[16 * 4];
LDPC_Encoder* LDPC_WiFi80211ad(const int rate, const int length, const bool lifts, int *seed 	)
{
	bool lift2 = lifts; // Turn on Samsung lifting for N=1344
	static const int lift12[16*8] = { 
		0, -1, 1, -1, 0, -1, 1, -1, 0, -1, -1, -1, -1, -1, -1, -1,
		0, -1, 0, -1, 1, -1, -1, 1, 0, 0, -1, -1, -1, -1, -1, -1,
		-1, 0, -1, 1, -1, 0, -1, 1, -1, 1, 0, -1, -1, -1, -1, -1,
		-1, 1, -1, 1, -1, 1, 0, -1, -1, -1, 0, 0, -1, -1, -1, -1,
		0, -1, 1, -1, 1, -1, 0, -1, 0, -1, -1, 1, 0, -1, -1, -1,
		1, -1, 1, -1, -1, 1, -1, 0, -1, 1, -1, 1, -1, 0, -1, -1,
		-1, 0, -1, 0, -1, 1, -1, 0, -1, -1, 0, -1, -1, 1, 0, -1,
		-1, 0, -1, 1, 0, -1, 0, -1, 0, -1, -1, -1, 1, -1, 0, 0 };

	static const int lift58[16*6] = {
		0,0,1,1,0,0,1,1,-1,1,0,-1,-1,-1,-1,-1,
		0,1,-1,1,-1,1,0,0,1,1,0,0,-1,-1,-1,-1,
		0,-1,1,-1,1,-1,0,-1,0,-1,-1,1,0,-1,-1,-1,
		1,-1,1,-1,-1,1,-1,0,-1,1,-1,1,1,0,-1,-1,
		-1,0,-1,0,-1,1,-1,0,-1,0,0,-1,-1,1,0,-1,
		-1,0,-1,1,0,-1,0,-1,0,-1,-1,-1,-1,-1,0,0 };	
	
	static const int lift34[16*4] = { 
		0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, -1, -1, -1,
		1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, -1, -1,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1,
		1, 0, 1, 1, 0, 1, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0 };



	static const int lift1316[16*3] = { 
		1,0,1,1,1,1,0,0,1,1,0,1,1,0,-1,-1,
		0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,-1,
		1,0,1,1,0,1,0,1,0,0,1,0,1,0,0,0 };


	static const int lift783[16 * 4] = {
		0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, -1, -1, -1, // Reg1
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1,	// Reg3
		1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, -1, -1,	// Reg 2
		1, 0, 1, 1, 0, 1, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0 };	// Reg 4
	//static const int lift783[16 * 4] = {
	//	0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, -1, -1, -1, // Reg1
	//	1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, -1, 1, 0, -1,	// Reg3
	//	1, 0, 1, 1, 1, 1, 0, 0, 1,  1, 0, 1, 1, 0, -1, -1,	// Reg 2
	//	0, 1, 0, 0, 0, 0, 1, 1, 0, -1, 1, 0, 0, 1, 0, 0 };	// Reg 4

	//	1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, -1, -1, -1, // Flip 1
	//	1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, -1, 0, 1, -1,  // Flip 3 
	//	0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, -1, -1,  // Flip 2
	//	0, 1, 0, 0, 1, 0, 1, 0, 1, -1, 0, 1, 0, 1, 1, 1 }; // Flip 4

	//static const int lift783[16 * 4] = {
	//	0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, -1, -1, -1,	//Reg 1
	//	1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, -1, 0, 1, -1,   // Flip 3
	//	1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, -1, -1,   // Reg 2
	//	0, 1, 0, 0, 1, 0, 1, 0, 1, -1, 0, 1, 0, 1, 1, 1 };  // Flip 4


	static const int lift784[16*4] = { 
		1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, -1, -1,
		0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1, 1, 0, -1,
		0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, -1, -1, -1,
		1, 0, 1, 1, 0, 1, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0 };

	const int *lift=0;


	int z;
	int *tab;
	int K;

	if (!lift2 && length % 672 != 0)
	{
		printf("LDPC_WiFi80211ad: the code length must be a multiple of 672,.\n");
		return 0;
	}
	else if (lift2 && (length != 672 && length != 1344))
	{
		printf("LDPC_WiFi80211ad: the code length must be  672 or 1344 for Samsung lifting.\n");
		return 0;
	}


	int N=length;


	int ncol;
	if (lift2 && N == 1344)
	{
		printf("LDPC_WiFi80211ad: Samsung lifting for generating length 1,344 code.\n");
		ncol = 32;
		z = 42;
	}
	else
	{
		ncol = 16;	
		z=42*(length/672);
		if(N>672)printf("LDPC_WiFi80211ad: Increasing z (%d) for generating length %d code.\n",z,N);

	}
	int ncheck = ncol;
	int nrow;
	int*sumrow = 0;
	switch(rate)
	{

	case 12:
		K=N/2;
		tab = WIFI80211ad_12EncTable;
		nrow = (N - K) / z;
		lift = lift12;
		break;

	case 58:
		K=(N*5)/8;
		tab =WIFI80211ad_58EncTable;
		nrow = (N - K) / z;
		lift = lift58;
		break;

	case 34:
		K = (N * 3) / 4;
		tab=WIFI80211ad_34EncTable;
		nrow = (N - K) / z;
		lift = lift34;

		break;

	case 1316:// 1316 or puncturing
	case 786:
		K = (N*13)/16;
		tab = WIFI80211ad_1316EncTable;
		nrow = (N - K) / z;
		lift = lift1316;
		break;

	case 781: // First two rows of 1316
		printf("Pick first two rows from 13/16 code generate rate 7/8 code.\n");
		K = (N*7)/8;
		tab = WIFI80211ad_781EncTable;
		nrow = (N - K) / z;
		break;

	case 782: // sum of first two rows from 1316
		printf("Sum pair of rows 12 from 13/16 code generate rate 7/8 code.\n");
		K = (N*7)/8;
		tab = WIFI80211ad_782EncTable;
		if (*seed == 0)
		{
			lift = lift1316;
		}
		else
		{
			int i;
			for (i = 0; i < 1 * 16; i++)
			{
				liftx[i] = unif_int(*seed) & 1;
			}
			for (; i < 2 * 16; i++)// This force the others to be complementary
			{
				liftx[i] = liftx[i - 16] ^ 1;
			}
			for (; i < 3 * 16; i++)
			{
				liftx[i] = unif_int(*seed) & 1;
			}
			lift = (const int*)liftx;
		}
		nrow = 3;
		sumrow = new int[nrow];
		sumrow[0] = 0;
		sumrow[1] = 0;
		sumrow[2] = 1;
		ncheck = 2 * ncol;
		break;

	case 783: // sum pairs of rows (13-24) from 34 (best)
		printf("Sum pairs of rows (13-24) from 34 code generate rate 7/8 code.\n");
		K = (N*7)/8;
		tab = WIFI80211ad_783EncTable;
		if (N == 1344)
		{
			if (*seed == 0)
			{
				lift = lift783;
			}
			else
			{
				int i;
	//			printf("Independent lifting.\n");
				printf("Identical lifting.\n");
	//			printf("Complementary lifting.\n");
				for (i = 0; i < 1* 16 ; i++)
				{
					liftx[i] = unif_int(*seed) & 1;
				}
				for (     ; i < 2*16 ; i++)// This force the others to be complementary
				{
	//				liftx[i] = unif_int(*seed) & 1;// This is random
					liftx[i] = liftx[i - 16] ^ 1;// This force the others to be complementary
	//				liftx[i] = liftx[i - 16];// This force the others to be the same
				}
				for (; i < 3*16; i++)
				{
					liftx[i] = unif_int(*seed) & 1;
				}
				for (; i < 4*16; i++)
				{
	//				liftx[i] = unif_int(*seed) & 1;// This is random
					liftx[i] = liftx[i - 16] ^ 1;// This force the others to be complementary
	//				liftx[i] = liftx[i - 16] ;// This force the others to be the same
				}

				lift = (const int*)liftx;
			}
		}
		nrow = 4;// Number of rows in the table
		sumrow = new int[nrow];
		sumrow[0] = 0;
		sumrow[1] = 0;
		sumrow[2] = 1;
		sumrow[3] = 1;
		ncheck = 2 * ncol;
		break;	

	case 784:// sum pairs of rows (23-14) from 34
		printf("Sum pairs of rows (23-14) from 34 code generate rate 7/8 code.\n");
		K = (N*7)/8;
		tab = WIFI80211ad_784EncTable;
		if (N == 1344)
		{
			if (*seed == 0)
			{
				lift = lift784;
			}
			else
			{
				int i;
				for (i = 0; i < 1 * 16; i++)
				{
					liftx[i] = unif_int(*seed) & 1;
				}
				for (; i < 2 * 16; i++)// This force the others to be complementary
				{
					liftx[i] = liftx[i - 16] ^ 1;
				}
				for (; i < 3 * 16; i++)
				{
					liftx[i] = unif_int(*seed) & 1;
				}
				for (; i < 4 * 16; i++)// This force the others to be complementary
				{
					liftx[i] = liftx[i - 16] ^ 1;
				}
				lift = (const int*)liftx;
			}
		}

		nrow = 4;
		sumrow = new int[nrow];
		sumrow[0] = 0;
		sumrow[1] = 0;
		sumrow[2] = 1;
		sumrow[3] = 1;
		ncheck = 2 * ncol;
		break;

	default:
		printf("LDPC_WiFi80211ad: Rate %d not available.\n",rate);
		exit(1);
		break;
	}

	int r=N-K;
	int* h=new int[(N-K)*(ncheck)];
	int*hh;
	int i,j,k,l,ii;
	int shift;
	if (!lift2 || N != 1344)
	{
		for(j=0;j<z;j++)
		{
			for(i=ii=0;i<r/z;i++)
			{
				l=0;
				hh = h + ncheck*(i*z + j);
				do
				{
					for(k=0;k<ncol;k++)
					{
						if(tab[ii*ncol+k]>=0)
						{
							shift=tab[ii*ncol+k];
							hh[l++] = k*z+ (shift+j)%z;
						}
					}
					ii++;
				} while (sumrow != 0 && sumrow[ii] == i);
				if(l<ncheck)hh[l]=-1;
			}

		}
	}
	else // Samsung lifting for N=1344
	{
		int k1, k2,lif,i1,i2;
		int dsh;
		dsh = 0;
//		dsh = z/2;
//		dsh = (unif_int(*seed) % z);
		for(j=0;j<z;j++)
		{
			for(i=ii=0;i<r/z;i++)// Loop on rows
			{
				l=0;
				hh = h + ncheck*(i*z + j);
				i1 = i / 2;
				i2 = i % 2;
				do
				{
					for(k=0;k<ncol;k++)
					{
						k1 = k / 2;
						k2 = k % 2;
						shift = tab[i1*(ncol / 2) + k1];
						lif   = lift[i1*(ncol / 2) + k1];
	//					printf("shift=%2d\tlift=%2d\n", shift, lif);
						if(shift>=0)
						{
							if ((i2^k2) == lif)hh[l++] = k*z + (shift + (k2^i1&1)*dsh + j) % z;
						}
					}
					i1++;
				} while (sumrow != 0 && sumrow[i1] == i/2);
				if(l<ncheck)hh[l]=-1;
			}

		}
	}
	delete[] sumrow;
	int x = 0;
	LDPC_Encoder* ldpc=new LDPC_Encoder;
	ldpc->x     = x;
	ldpc->z     = z;
	ldpc->K     = K;
	ldpc->N	    = ldpc->N1    =N;
	ldpc->Hconn = h;
	ldpc->ncheck = ncheck;
	if (rate == 786)
	{
		printf("Applying puncturing to generate rate 7/8 code.\n");
		int* pattern = new int[N];

		for (i = 0; i < N; i++)pattern[i] = 1;
		for (i = K; i < K + 48*(N/672); i++)pattern[i] = 0;
		ldpc->SetPuncturing(N,pattern);
	}
	ldpc->SetEncodingMethod(0);
	return ldpc;
}
#include "cppincludes/WiMAXCodes.cpp"

LDPC_Encoder* LDPC_WIMAX(const int rate, const int N, const bool A)
{
	if(N%24!=0)return 0;
	int z=N/24;
	if(z<24 ||z%4!=0)
	{
		printf("LDPC_WIMAX: code length %d not allowed.\n",N);
		exit(1);
	}
	int *tab;
	int r,x;
	int ncheck=24;
	int K;
	switch(rate)
	{
	case 12:
		tab=WIMAX_12EncTable;
		r=12;

		K=N/2;
		x  =0;
		break;
	case 23:
		if(A)
		{
			tab=WIMAX_23AEncTable;
			x  =0;
		}
		else 
		{
			tab=WIMAX_23BEncTable;
			x  =0;
		}
		r=8;
		K=(N*2)/3;
		break;
	case 34:
		if(A)
		{
			tab = WIMAX_34AEncTable;
			x   = 0;
		}
		else 
		{
			tab= WIMAX_34BEncTable;
			x  = (80*z)/96;
		}
		r=6;
		K=(N*3)/4;
		break;
	case 56:
		tab=WIMAX_56EncTable;
		r=4;
		K=(N*5)/6;
		x  = 0;
		break;
	default:
		printf("LDPC_WIMAX: Rate %d not allowed.\n",rate);
		exit(1);

	}

	int* h=new int[(N-K)*ncheck];
	int*hh;
	int i,j,k,l;
	int shift;
	l=0;
	for(i=0;i<r;i++)
	{
		for(j=0;j<z;j++)
		{
			hh=h + 24*(i*z+j);
			l=0;
			for(k=0;k<24;k++)
			{
				if(tab[i*24+k]>=0)
				{
					if(rate==23 && A)
					{
						shift = tab[i*24+k] % z;
					}
					else
					{
						shift=(int)(((double)tab[i*24+k]*z)/96.);
					}
					hh[l++] = k*z+ (shift+j)%z;
				}

			}
			hh[l]=-1;
		}
		
	}

	LDPC_Encoder* ldpc=new LDPC_Encoder;
	ldpc->x     = x;
	ldpc->z     = z;
	ldpc->K     = K;
	ldpc->N	    = ldpc->N1    =N;
	ldpc->Hconn = h;
	ldpc->ncheck = ncheck;
	ldpc->SetEncodingMethod(2);
	return ldpc;

}
int* Gaussian_Elimination(unsigned int *H, const int n, const int k)
{
	int i,j,l;
	unsigned int* h;
	unsigned int* h2;
	int* p=new int[n];
	for(i=0;i<n;i++)p[i]=i;
	int temp;
	int M=n/32;	//elements in a row
	if(n%32!=0)M++;

	for(i=0;i<k;i++)// 
	{
		// Find a column j with a 1 in i-th position
		h=H+ i*M;
		j=i;
		while(getbit(h,p[j])==0){j++;}

		// Swap i and j columns
		temp=p[i];
		p[i]=p[j];
		p[j]=temp;

		// Zeroes following rows
		for(j=i+1;j<k;j++)
		{
			h2=H+j*M;
			if(getbit(h2,p[i]))for(l=0;l<M;l++)h2[l]^=h[l];
		}

		// Zeroes previous rows
		for(j=i-1;j>=0;j--)
		{
			h2=H+j*M;
			if(getbit(h2,p[i]))for(l=0;l<M;l++)h2[l]^=h[l];
		}
	}
	return p;
}
int* LDPC_Encoder::PEG(const int K, const int N, const int nedges, const int *degv, const int* degc)
{
	const int lmax = 10000;
	int i, j;
	int *vdeg = new int[N];
	int *cdeg = new int[N - K];
	int *vcon = new int[nedges];
	int *ccon = new int[nedges];
	int **vvcon = new int*[N];
	int **cccon = new int*[N - K];

	vvcon[0] = vcon;
	for (i = 1; i < N; i++)vvcon[i] = vvcon[i - 1] + degv[i - 1];
	cccon[0] = ccon;
	for (i = 1; i < N - K; i++)cccon[i] = cccon[i - 1] + degc[i - 1];

	int *I = new int[N-K];

	for (i = 0; i < N; i++)    vdeg[i] = 0;
	for (i = 0; i < N - K; i++)cdeg[i] = 0;
	int dmin, cbest;
	int c, v;
	bool again;

	int l, k1, k2;
	for (i = 0; i < N; i++)// Loop on VN
	{
		for (j = 0; j < N-K; j++)I[j] = lmax;

		fprintf(stdout, "VN %5d (%d):\n", i,degv[i]);
		for (j = 0; j < degv[i]; j++)// Loop on VN node degree
		{
			if (j > 0)
			{
				again = true;
				for (l = 1; again; l++)// Subgraph expansion
				{
					again = false;
					// Check node expansion
					for (c = 0; c < N - K; c++)
					{
						if (I[c] == lmax)continue;
						for (k1 = 0; k1 < cdeg[c]; k1++)
						{
							v = cccon[c][k1];  // Connected variable node
							for (k2 = 0; k2<vdeg[v]; k2++)
							{
								if (I[vvcon[v][k2]]>l)
								{
									I[vvcon[v][k2]] = l;
									again = true;
								}
							}

						}
					}
				}
			}

			// Search node with highest girth and smallest degree
			dmin = 1000;
			l = 0;
			for (c = 0; c < N-K; c++)
			{
				if (I[c] > l)
				{
					l = I[c];
					dmin = cdeg[c];
					cbest = c;
				}
				else if (I[c] == l)
				{
					if (cdeg[c] < dmin)
					{
						dmin = cdeg[c];
						cbest = c;
					}
				}
			}

			// Look for the first with largest girth and smallest dmin
			c = cbest;
			vvcon[i][vdeg[i]] = c;
			vdeg[i]++;
			cccon[c][cdeg[c]] = i;
			cdeg[c]++;
			I[c] = 0;
			fprintf(stdout, "\t\t%4d", c);
			fprintf(stdout, "->%d %d\n", dmin, l);

		}
	}
	return 0;
}
void LDPC_Encoder::AddDegreeOneVN(
	const int Ndegone
	)
{
	int i,j,row,col,temp;
	int M = N1 - K1;
	int* Hconnnew = new int[ncheck*(M+ Ndegone)];
	for (i = 0; i < ncheck*M; i++)Hconnnew[i] = Hconn[i];
	delete[] Hconn;
	Hconn = Hconnnew;

	int ncheck2 = ncheck;
//	printf("%d\n",ncheck2);

	int *H1       = Hconn+ncheck*M;
	int l = -1;
	for (i = 0; i < Ndegone; i++)
	{
		for (j = 0; j < ncheck2-1; j++)
		{
			//l+=23;
			//H1[i*ncheck + j] = l%K1;
			do {
				l++;
				row = l%M;
				col = (l / M) % ncheck;
				temp = Hconn[row*ncheck + col];
			} while (temp < 0);
			H1[i*ncheck + j] = temp;
		}
		H1[i*ncheck + ncheck2 - 1] = N1 + i;
		if(ncheck2<ncheck)H1[i*ncheck + ncheck2] = -1;
	}
	N1 += Ndegone;
	N  += Ndegone;
}

//int  SetMultiEdge_old(
//	const int N,
//	const double* NU,
//	const double* MU
//	)
//{
//	int success = 1;
//	int i, j, k, ipick, seed, temp;
//	double R = EvalMN(NU, 0, 0) - EvalMN(MU, 0, 0);
//	this->N = N1 = N;
//	this->K = K1 = (int)(N*R+0.5);
//	int ntypes = (int)NU[1];
//	int *perm1 = 0, *perm2 = 0, *perm3 = 0, *perm4 = 0;
//	int *imaxN = 0, *counter1 = 0, *counter2 = 0, *base = 0,*ldegn=0,*rdegn=0;
//	int *nones1 = new int[ntypes];
//	int *nones2 = new int[ntypes];
//	for (j = 0; j < ntypes; j++)nones1[j] = 0;
//	for (j = 0; j < ntypes; j++)nones2[j] = 0;
//	int l, l2;
//	int m = 0;
//	double cum = 0.;
//	int d;
//	int* NN = new int[ntypes]; // Number of Nodes connected to edge i
//	int* MM = new int[ntypes]; // Number of Check nodes connected to edge type i
//	int** degseq = new int*[ntypes];
//	for (i = 0; i < ntypes; i++)
//	{
//		degseq[i] = new int[N];
//		NN[i] = 0;
//		MM[i] = 0;
//	}
//
//	int ndv = (int)NU[0];
//	ldegn = new int[ndv]; // Number of nodes of type i
//	m = 0;
//	cum = 0.;
//	for (i = 0;i<ndv - 1; i++)
//	{
//		cum += NU[2+i*(ntypes+1)];
//		d = (int)(cum*N + 0.5) - m;
//		m += ldegn[i] = d;
//	}
//	ldegn[ndv - 1] = N - m;
//	for (i = 0; i < ndv; i++)
//	{
//		for (j = 0; j < ntypes; j++)
//			nones1[j] += (int)NU[2+i*(ntypes+1)+1+j]*ldegn[i];
//	}
//
//	printf("Finite size Multiedge (%d) degree  distribution\n", N);
//	for (i = 0; i < ndv; i++)
//	{
//		fprintf(stdout, "NU[");
//		for (j = 0; j < ntypes - 1; j++)
//		{
//			fprintf(stdout,"%2d,", (int)NU[2 + i*(ntypes + 1) + 1 + j]);
//			if (NU[2 + i*(ntypes + 1) + 1 + j]>0)NN[j] += ldegn[i];
//		}
//		fprintf(stdout, "%2d]=", (int)NU[2 + i*(ntypes + 1) + 1 + j]);
//		if (NU[2 + i*(ntypes + 1) + 1 + j]>0)NN[j] += ldegn[i];
//		fprintf(stdout, "%d\n", ldegn[i]);
//	}
//	for (j = 0; j < ntypes; j++)
//		fprintf(stdout, "NN[%d]=%d\n", j, NN[j]);
//
//	int ndc = (int)MU[0];
//	rdegn = new int[ndc]; // Number of nodes of type i
//	m = 0;
//	cum = 0.;
//	for (i = 0; i<ndc - 1; i++)
//	{
//		cum += MU[2 + i*(ntypes + 1)];
//		d = (int)(cum*N + 0.5) - m;
//		m += rdegn[i] = d;	
//	}
//	rdegn[ndc - 1] = N - K - m;
//	int deg;
//	int maxcheck;
//	maxcheck = 0;
//	for (i = 0; i < ndc; i++)
//	{
//		fprintf(stdout, "MU[");
//		deg = 0;
//		for (j = 0; j < ntypes - 1; j++)
//		{
//			deg += MU[2 + i*(ntypes + 1) + 1 + j];
//			fprintf(stdout, "%2d,", (int)MU[2 + i*(ntypes + 1) + 1 + j]);
//			if (MU[2 + i*(ntypes + 1) + 1 + j]>0)MM[j] += rdegn[i];
//		}
//		if (rdegn[i] > 0 && deg>maxcheck)maxcheck = deg;
//		fprintf(stdout, "%2d]=", (int)MU[2 + i*(ntypes + 1) + 1 + j]);
//		if (MU[2 + i*(ntypes + 1) + 1 + j]>0)MM[j] += rdegn[i];
//		fprintf(stdout, "%d (deg=%d)\n", rdegn[i],deg);
//	}
//	for (j = 0; j < ntypes; j++)
//		fprintf(stdout, "MM[%d]=%d\n", j,MM[j]);
//
//
//	fprintf(stdout, "Maximum check degree= %d\n", maxcheck);
//	maxcheck+=2;
//	for (i = 0; i < ndc; i++)
//	{
//		for (j = 0; j < ntypes; j++)
//			nones2[j] += (int)MU[2 + i*(ntypes + 1) + 1 + j] * rdegn[i];
//	}
//
//
//	int nontot;
//	counter1 = new int[ntypes];
//	counter2 = new int[ntypes];
//	base    = new int[ntypes];
//	for (j = nontot= 0; j < ntypes; j++)
//	{
//		base[j] = nontot;
//		fprintf(stdout, "%d\t%d\t%d\n", nones1[j], nones2[j], nones1[j] - nones2[j]);
//		nontot += (nones1[j]+nones2[j])/2; // The number of ones is the average
//	}
//
//	perm1 = new int[nontot]; // Perm grouping edge types at variable node 
//	perm2 = new int[nontot]; // Perm grouping edge types at check nodes  
//	perm3 = new int[nontot]; // Perm of edge type  
//	perm4 = new int[nontot]; // Perm of edge type  
//
//	int imaxM;
//	int a,ll;
//
//	// Construct permutation (perm1) for storing sockets in order of edge type
//	imaxN = new int[ntypes];
//
//	int D;
//	for (j = 0; j < ntypes; j++)// Imax is the node with maximum variable degree of that type
//	{
//		NN[j] = 0;
//		counter1[j] = 0;
//		if (nones1[j] < nones2[j])D = abs(nones1[j] - nones2[j]) / 2;
//		else                      D = (int)(fabs(nones1[j] - nones2[j]) / 2. + 0.9999);
//		if (D == 0)continue;
//		imaxN[j] = -1;
//		for (k = 0; k < ndv; k++)// Compute the maximum degree for that type
//		{
//			if (ldegn[k] >= D && NU[2 + k*(ntypes + 1) + j + 1]>0)
//			{
//				if(imaxN[j]==-1)imaxN[j] = k;
//				else           
//					if(NU[2 + k*(ntypes + 1) + j + 1]>NU[2 + imaxN[j]*(ntypes + 1) + j + 1])imaxN[j] = k;
//			}
//		}
//		if (imaxN[j] ==-1)
//		{
//			printf("Failure in ME LDPC code construction (%d).\n",j);
//			success = 0;
//			goto exit;
//		}
//	}
//	for (k = ll = 0; k < ndv; k++)// Loop on node types
//	{
//		for (i = 0; i < ldegn[k]; i++)// Number of nodes of given degree types
//		{
//			for (j = 0; j < ntypes; j++)
//			{
//				if (NU[2 + k*(ntypes + 1) + j + 1] == 0)continue;
//				a = 0;
//				if (k == imaxN[j])
//				{
//					if(nones1[j] < nones2[j])D = abs(nones1[j] - nones2[j]) / 2;
//					else                     D = (int)(fabs(nones1[j] - nones2[j]) / 2. + 0.9999);
//					if (i >= ldegn[k] - D)
//					{
//						if (nones1[j] > nones2[j])a = -1;
//						else                      a = 1;
//					}
//				}
//				for (l = 0; l < NU[2 + k*(ntypes + 1) + j + 1]+a; l++)
//				{
//					perm1[base[j] + counter1[j]] = ll;
//					counter1[j]++;
////					ll++;  // This is for storing the edge index
//				}
//				degseq[j][NN[j]++] = NU[2 + k*(ntypes + 1) + j + 1] + a;
//			}
//			ll++; // This is for storing node index
//		}
//	}
//
//
//	// PEG permutations on each type with almost regular right degree
//
//	BigGirth *big = 0;
//	for (j = 0; j < ntypes; j++)
//	{
//		big = new BigGirth(MM[j],NN[j],degseq[j],0,0);
//		big->
//
//	}
//
//	// Construct permutation for storing sockets at check nodes
//	for (j = 0; j < ntypes; j++)// Imax is the node with variable degree of that type
//	{
//		counter2[j] = 0;
//		if (nones1[j] > nones2[j])D = abs(nones1[j] - nones2[j]) / 2;
//		else                      D = (int)(fabs(nones1[j] - nones2[j]) / 2. + 0.9999);
//		if (D == 0)continue;
//		imaxN[j] = -1;
//		for (k =0; k < ndc; k++)// Compute the maximum degree for that type
//		{
//			if (rdegn[k] >= D && MU[2 + k*(ntypes + 1) + j + 1]>0)
//			{
//				if(imaxN[j]==-1)imaxN[j] = k;
//				else
//					if (MU[2 + k*(ntypes + 1) + j + 1]>MU[2 + imaxN[j] * (ntypes + 1) + j + 1])imaxN[j] = k;
//			}
//		}
//		if (imaxN[j] == -1)
//		{
//			printf("Failure in code construction (%d).\n",j);
//			success = 0;
//			goto exit;
//		}
//	}
//	delete[] Hconn;
//	Hconn = new int[maxcheck*(N - K)];
//	int* H = Hconn;
//	int vnode;
//	for (k = ll= 0; k < ndc; k++)// Loop on degree types
//	{
//		for (i = 0; i < rdegn[k]; i++)// Number of nodes of given degree types
//		{
//			m = 0;
//			for (j = 0; j < ntypes; j++)// Loop on edge types
//			{
//				if (MU[2 + k*(ntypes + 1) + j + 1] == 0)continue;
//				a = 0;
//				if (k == imaxN[j])
//				{
//					if (nones1[j] > nones2[j])D = abs(nones1[j] - nones2[j]) / 2;
//					else                      D = (int)(fabs(nones1[j] - nones2[j]) / 2. + 0.9999);
//
//					if (i >= rdegn[k] - D)
//					{
//						if (nones1[j]>nones2[j])a =  1;
//						else                    a = -1;
//					}
//				}
//
//				for (l = 0; l < MU[2 + k*(ntypes + 1) + j + 1]; l++)// Loop on sockets of given types
//				{
//					perm3[ll] = base[j] + counter2[j];
//					vnode = perm1[perm2[perm3[ll]]];
//					H[m++] = vnode;
//					ll++;
//					counter2[j]++;
//				}
//			}
//			if (m < maxcheck)H[m] = -1;
//			H += maxcheck;
//		}
//	}
//
//	for (j = 0; j < ntypes; j++)
//	{
//		printf("counter1[%d]=%d,counter2[%d]=%d\n",j,counter1[j],j,counter2[j]);
//		if (counter1[j] != counter2[j])
//		{
//			printf("Error in allocation of sockets on type %d.\n",j);
//			delete[] Hconn;
//			success = 0;
//			goto exit;
//		}
//	}
//
//
//	///* Combine permutations (Readind order for check nodes sockets )*/
//	for (j = 0; j < nontot; j++)
//	{
//		perm4[j] = perm1[perm2[perm3[j]]];
//	}
//
//	// Now write the parity check matrix
//exit:
//	if (success)
//	{
//		this->K = K;
//		this->K1 = K1;
//		this->N = N;
//		this->N1 = this->N = N;
//		this->ncheck = maxcheck;
//		this->shorten = 0;
//		this->gap = 0;
//	}
//	delete[] nones1;
//	delete[] nones2;
//	delete[] perm1;
//	delete[] perm2;
//	delete[] perm3;
//	delete[] perm4;
//	delete[] imaxN;
//	delete[] counter1;
//	delete[] counter2;
//	delete[] base;
//	delete[] ldegn;
//	delete[] rdegn;
//	return success;
//}

//int  LDPC_Encoder::SetMultiEdge(
//	const int N,
//	const double* NU, // Left multi degree multi profile
//	const double* MU, // Right quasi regular profile (Only binary types are allowed)
//	const int Z
//	)
//{
//	int success = 1;
//	int i, j, k;
//	double R = EvalMN(NU, 0, 0) - EvalMN(MU, 0, 0);
//	this->N = N1 = N;
//	this->K = K1 = (int)(N*R + 0.5);
//	int ntypes = (int)NU[1];
//	int *ldegn = 0, *rdegn = 0;
//	int *nones = new int[ntypes];
//	for (j = 0; j < ntypes; j++)nones[j] = 0;  // Number of edges of each type
//	int m = 0;
//	double cum = 0.;
//	int d;
//
//	int* degseq = new int[N];
//
//	int ndv = (int)NU[0];
//	ldegn = new int[ndv]; // Number of nodes of type i
//	m = 0;
//	cum = 0.;
//	for (i = 0; i < ndv - 1; i++)
//	{
//		cum += NU[2 + i*(ntypes + 1)];
//		d = Z*(int)(cum*(N/Z) + 0.5) - m;
//		m += ldegn[i] = d;
//	}
//	cum += NU[2 + i*(ntypes + 1)];
//	ldegn[ndv - 1] = N - m;
//	for (i = 0; i < ndv; i++)
//	{
//		for (j = 0; j < ntypes; j++)
//			nones[j] += (int)NU[2 + i*(ntypes + 1) + 1 + j] * ldegn[i];
//	}
//
//	printf("Finite size Multiedge (%d) left degree  distribution\n", N);
//	for (i = 0; i < ndv; i++)
//	{
//		fprintf(stdout, "NU[");
//		for (j = 0; j < ntypes - 1; j++)
//		{
//			fprintf(stdout, "%2d,", (int)NU[2 + i*(ntypes + 1) + 1 + j]);
//		}
//		fprintf(stdout, "%2d]=", (int)NU[2 + i*(ntypes + 1) + 1 + j]);
//		fprintf(stdout, "%d\n", ldegn[i]);
//	}
//
//	int* MM = new int[ntypes]; // Number of Check nodes connected to edge type i
//	for (i = 0; i < ntypes; i++)MM[i] = 0;
//	int ndc = (int)MU[0];
//	rdegn = new int[ndc]; // Number of nodes of type i (only 1 or zero is possible)
//	m = 0;
//	cum = 0.;
//	for (i = 0; i<ndc - 1; i++)
//	{
//		cum += MU[2 + i*(ntypes + 1)];
//		d = Z*(int)(cum*(N/Z) + 0.5) - m;
//		m += rdegn[i] = d;
//		for (j = 0; j < ntypes; j++)
//		{
//			if (MU[2 + i*(ntypes + 1) + 1 + j]>0)MM[j] += rdegn[i];
//		}
//	}
//	cum += MU[2 + i*(ntypes + 1)];
//	rdegn[ndc-1] = N - K - m;
//	for (j = 0; j < ntypes; j++)
//	{
//		if (MU[2 + (ndc-1)*(ntypes + 1) + 1 + j]>0)MM[j] += rdegn[ndc-1];
//	}
//
//
//
//	for (i = 0; i < ndc; i++)
//	{
//		fprintf(stdout, "MU[");
//		for (j = 0; j < ntypes - 1; j++)
//		{
//			fprintf(stdout, "%2d,", (int)MU[2 + i*(ntypes + 1) + 1 + j]);
//		}
//		fprintf(stdout, "%2d]=", (int)MU[2 + i*(ntypes + 1) + 1 + j]);
//		fprintf(stdout, "%d \n", rdegn[i]);
//	}
//	for (j = 0; j < ntypes; j++)
//		fprintf(stdout, "MM[%d]=%d\n", j, MM[j]);
//
//	
//	// Construct permutations
//	int maxcheck;
//
//	BigGirth *big = 0;
//	int ll;
//	int** Hc = new int*[ntypes];
//	int* mc = new int[ntypes];
//	maxcheck = 0;
//	for (j = 0; j <ntypes; j++)
//	{
//		for (k = ll = 0; k < ndv; k++)// Loop on node types
//		{
//			for (i = 0; i < ldegn[k]/Z; i++)// Number of nodes of given degree types
//			{
//				degseq[ll++] = (int)NU[2 + k*(ntypes + 1) + j + 1];
//			}
//		}
//		big = new BigGirth(MM[j],N, degseq,1, 100,Z,true,0,false,127560);
////		big->writeToFile_Hcompressed("Hcomp.txt", Z);
////		printf("\n Completed.\n");
//		//		delete[] Gconn; Gconn = 0;
//		mc[j] = 0;
//		Hc[j] = big->Build_Hcompressed(mc[j]);
//		maxcheck += mc[j];
//		delete big;
//	}
//
//	Hconn = new int[(N-K)*maxcheck];
//	int* p = new int[ntypes];
//	for (j = 0; j < ntypes; j++)p[j] = 0;
//	int m1 ;
//
//	// Merging  check nodes 
//	for (k = ll = 0; k < ndc; k++)// Loop on check node types (binary (1,1,0)...)
//	{
//		for (i = 0; i < rdegn[k]; i++)// Number of nodes of given degree type (merging multiple types)
//		{
//			m1 = 0;
//			for (j = 0; j < ntypes; j++)
//			{
//				if (MU[2 + k*(ntypes + 1) + j + 1] == 0)continue;
//				for (m = 0; m < mc[j]; m++)
//				{
//					if (Hc[j][p[j] * mc[j] + m] == -1) break;
//					Hconn[ll*maxcheck + m1] = Hc[j][p[j]*mc[j]+m];
//					m1++;
//				}
//				p[j]++;
//				
//			}
//			if (m1 < maxcheck)Hconn[ll*maxcheck + m1] = -1;
//			ll++; // This is for storing node index
//		}
//	}
//
//
//
//	// Now write the parity check matrix
//	if (success)
//	{
//		this->K = K;
//		this->K1 = K1;
//		this->N = N;
//		this->N1 = this->N = N;
//		this->ncheck = maxcheck;
//		this->shorten = 0;
//		this->gap = 0;
//	}
//	delete[] nones;
//	delete[] ldegn;
//	delete[] rdegn;
//	return success;
//}