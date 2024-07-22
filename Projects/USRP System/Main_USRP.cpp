#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <mutex>
#include "tool_functions.h"
#include "Parameter_Manager.h"
#include "PN_Source.h"
#include "LDPC_decoder.h"
#include "Modulator.h"
#include "Filter.h"
#include "AWGN_Channel.h"
#include "Time_Synchronizer.h"
#include "Demodulator.h"
#include "BER_meter.h"
#include "Spectrum.h"
#include "Speed.h"
#include "Automatic_Gain_Control.h"
#include "Pilot.h"
#include "Equalizer.h"
#include "Frame_sync.h"
#include "Tuning_centroids.h"
#include "USRP_Manager.h"

double RXpower;
double gainT = 15.;
double gainR = 28;
const int ff = 6;
const int sfifo=1<<ff;
inline int cd(const int a, const int b) 
{
	const int N = 1<<ff;
	int dm=((int)(b - a))%N;
	if(dm<0)dm+=N;
	return dm;
}
double* TX;
double* RX;
USRP_Manager* USRP = 0;
int pt1;
int pt2;
int pr1;
int pr2;
int nofsamp;
mutex mtx;

void SyncTX(__int64 &avewaittx, int &losttx)
{
	int dxt;
	__int64 wait;
	dxt = cd(pt2, pt1);
	if (dxt == sfifo - 1)
	{
		printf("PC : Lost TX packet (%d,%d)\n",  pt1, pt2);
		losttx++;
	}
	wait = 0;
	while (dxt >sfifo - 3)
	{
		mtx.lock();
		dxt = cd(pt2, pt1);
		mtx.unlock();
		wait++;
	}   // Avoid filling the same buffer
	avewaittx += wait;
}
void SyncRX(__int64 &avewaitrx, int &lostrx)
{
	int dxr;
	__int64 wait;
	dxr = cd(pr1, pr2);
	if (dxr == 0)
	{
		printf("PC : Lost RX packet (%d,%d)\n", pr1, pr2);
		lostrx++;
	}
	wait = 0;
	while (dxr <2)// Synch with USRP buffer
	{
		mtx.lock();
		dxr = cd(pr1, pr2);
		mtx.unlock();
		wait++;
	}
	avewaitrx += wait;
}
void Send()
{
	int res;
	while (pt2>=0)
	{
		res = USRP->TX(TX + (pt2 &(sfifo - 1)) * 2 * nofsamp);
		mtx.lock();
		pt2++;
		mtx.unlock();
	}
}
void Receive()
{
	int res;
	double *p;
	while (pr2>=0)
	{
		p = RX + (pr2&(sfifo - 1)) * 2 * nofsamp;
		res=USRP->RX(p);
		mtx.lock();
		pr2++;
		mtx.unlock();
	}
}
void TXRX()
{
	int res;
	double*p;
	while (pt2>=0)
	{
		p = TX + (pt2&(sfifo - 1)) * 2 * nofsamp;
		res=USRP->TX(p);

		p = RX + (pr2&(sfifo - 1)) * 2 * nofsamp;
		res=USRP->RX(p);

//		AGC(p);

		mtx.lock();
		pt2++;
		pr2++;
		mtx.unlock();
	}
}
int main()
{
	int N,D, H, P;
	int fmax;
	char nameinput[128]="Data/input_USRP.txt";
	double alpha;
	double decw;         // If different it alters the power computation
	double alphap;       // 0.005;
	double decwp;        //0.8;
	double afoff;
	double alphact;
	int coh = 5;
	bool specmeas = false;// Turn this on to turn on spectrum measures
	int Nd = 3;
	int Nfil = 10;
	int L0;
	double rolloff;
	int m,ns,nsymb,npsymb,nsx;
	double carrierFrequency;
	int actTX, actRX;
	double Rs,dfnorm,gamma;

	int ncand = 10;
	char name[132];
	double nu;
	TX = 0;
	RX = 0;

	Parameter_Manager a(40);
	FILE *input;
	input=fopen(nameinput,"r+");
	if(input==NULL)input=fopen(nameinput,"w+");
	a.Add_Parameter_From_File(input, "Name of output file", name);
	a.Add_Parameter_From_File(input, "Upsampling factor", &ns);
	a.Add_Parameter_From_File(input, "Number of frames", &fmax);

	a.Add_Parameter_From_File(input, "Header length (H)", &H);
	a.Add_Parameter_From_File(input, "Number of data sections in a frame (N)", &N);
	a.Add_Parameter_From_File(input, "Data Length (D)", &D);
	a.Add_Parameter_From_File(input, "Pilot length (P)", &P);

	a.Add_Parameter_From_File(input, "Maximum frequency offset [nu/Rs]", &nu);
	a.Add_Parameter_From_File(input, "SRRC filter roll-off", &rolloff);
	a.Add_Parameter_From_File(input, "QAM modulation efficiency", &m);
	a.Add_Parameter_From_File(input, "Upsampling factor", &ns);
	a.Add_Parameter_From_File(input, "AGC bandwidth", &gamma);
	a.Add_Parameter_From_File(input, "Observation window of timing", &L0);
	a.Add_Parameter_From_File(input, "Updating step of frequency estimator", &afoff);
	a.Add_Parameter_From_File(input, "Updating step of Equalizer", &alpha);
	a.Add_Parameter_From_File(input, "Updating step of PLL", &alphap);
	a.Add_Parameter_From_File(input, "Weight of decision aided error for Equalizer", &decw);
	a.Add_Parameter_From_File(input, "Weight of decision aided error for", &decwp);
	a.Add_Parameter_From_File(input, "Updating step of CT", &alphact);
	//USRP 
	a.Add_Parameter_From_File(input, "Activate TX?", &actTX);
	a.Add_Parameter_From_File(input, "Activate RX?", &actRX);
	a.Add_Parameter_From_File(input, "Spacing TX/RX carriers", &dfnorm);
	a.Add_Parameter_From_File(input, "USRP: Carrier frequency", &carrierFrequency);
	a.Add_Parameter_From_File(input, "USRP: Symbol rate", &Rs);
	a.Add_Parameter_From_File(input, "USRP: TX gain [dB]", &gainT);
	a.Add_Parameter_From_File(input, "USRP: RX gain [dB]", &gainR);
	

	fclose(input);
	if(a.error)
	{
		printf("Main: The file %s was not ready. Fill up all the required parameters.\n",nameinput);
		exit(1);
	}

	FILE* file = fopen(name, "w");
	FILE* out1 = fopen("Data/IQplot.txt", "w");
	FILE* out2 = fopen("Data/scatdg_BEFORE.txt", "w");
	FILE* out3 = fopen("Data/scatdg_MMSE.txt", "w");
	FILE* out4 = fopen("Data/Centroids.txt", "w");

	a.Print_Parameters_on_Line(file, true);
	fprintf(file, "gain TX\t");
	fprintf(file, "gain RX\t");
	fprintf(file, "Lost PC TX\t");
	fprintf(file, "Lost PC RX\t");
	fprintf(file, "Lost USB TX\t");
	fprintf(file, "Lost USB RX\t");
	fprintf(file, "\n");
	fflush(file);
	//	Open USRP and set carrrier frequency and IQ rate and buffer size

	do{
		if (actTX == 0 && actRX == 0)continue;
		if (Nd % 2 != 1)continue;
		if (nu == 0.)coh = H / 2;
		else         coh = (int)(1. / fabs(10.*nu)+0.5);
		if (coh > H / 2)coh = H / 2;
		if (coh < 1)coh = 1;

		a.Print_Parameters();
		nsymb  = N*D;
		npsymb = H + N*(D + P);
		nofsamp = ns*npsymb;

		// Check USRP presence and configure it	
		USRP = new USRP_Manager(actTX, actRX, Rs*ns, carrierFrequency, nofsamp, gainT, gainR);
		if (USRP->numberOfDevicesInSystem == 0 || USRP->nerr > 0)
		{
			printf("NIUSRP device is not available.\n");
			system("Pause");
			return 1;
		}
		if (dfnorm >0)// Set RX frequency to
		{
			printf("Set TX freq to %e\n", carrierFrequency + dfnorm*Rs);
			USRP->SetCarrier(1, carrierFrequency + dfnorm*Rs);
		}	
		else if (dfnorm<0)// Set RX frequency to
		{
			printf("Set RX freq to %e\n", carrierFrequency - dfnorm*Rs);
			USRP->SetCarrier(0, carrierFrequency - dfnorm*Rs);
		}


		Speed speed;
		PN_Source *Source;

		Modulator *MOD = 0;
		if(m==3)MOD= PSK_Modulator(m);
		else    MOD= QAM_Modulator(m);


		Pilot* PILOT = new Pilot;
		PILOT->SetRegularFrame(N, H, D, P); // Frame structure with header of length H and N repetitions of D+P
		int    *dd    = new int[m*(H + N*P)]; //source bits
		double* pilot = new double[2 * (H + N*P)];
		Source = new PN_Source((int)(log((H + N*P)*m) / log(2.) + 0.9999));
		Source->Run((H + N*P)*m, dd);
		MOD   ->Run((H + N*P), dd, pilot);
		PILOT->SetPilot((H + N*P), pilot);
		delete[] dd;
		delete[] pilot;
		delete Source;

		Source = new PN_Source;

		Filter* TXfil = SRRC(rolloff,ns,Nfil,true);
		TXfil->Set_Unitary_Energy();

		Filter* RXfil = SRRC(rolloff, ns, Nfil, true);
		RXfil->Set_Unitary_Energy();
		RXfil->SetMatched();

		// AGC
		Automatic_Gain_Control* AGC = new Automatic_Gain_Control;
		AGC->SetParameters(1, gamma);
		// Timing
		Time_Synchronizer* TSYNC = new Time_Synchronizer;
		TSYNC->SetParameters(ns,L0);

		// Frame Sync and frequency offset estimation
		Frame_Sync* FSYNC_FOFF = new Frame_Sync;
		double* ppp = PILOT->GetPilotSequence(H+N*(D+P));
		FSYNC_FOFF->SetParameters(npsymb,H, ppp, ncand,coh,afoff);
		FSYNC_FOFF->pointer+= Nd / 2;

		// MMSE_PLL 
		Equalizer* MMSE_PLL = new Equalizer;
		MMSE_PLL->SetParameters(MOD, Nd, alpha,1,0,0.,decw);
		MMSE_PLL->SetTraining(npsymb,ppp);
		MMSE_PLL->ActivatePLL(alphap, decwp);
		delete[] ppp;

		// Centroid tuning 
		Tuning_Centroids* CT = new Tuning_Centroids;
		CT->SetParameters(MOD,alphact);

		Demodulator* DEM = new Demodulator;
		Modulator*RXMOD = new Modulator;
		RXMOD->SetParameters(MOD->m, MOD->constellation, MOD->mapping);
		DEM->SetParameters(RXMOD);
		double SNR = 5;
		DEM->SetSigma(pow(10., -SNR / 10) / 2.);

		Delay* DEL = new Delay;
		DEL->SetParameters(N*m*D, 0);

		BER_meter* BER = new BER_meter[1];
		BER[0].SetParameters(0, 30);
		BER[0].SetFrameSize(N*D*m);
		BER[0].SetSoft();


		Spectrum* SP = new Spectrum[3];
		SP[0].SetParameters(10, 1.  / ns, "Data/spectrumTX.txt");
		SP[1].SetParameters(10, 1.  / ns, "Data/spectrumRX.txt");
		SP[2].SetParameters(10, 1. / ns, "Data/spectrumRX_afterSRRC.txt");


		/* Allocate buffers for storing input and outputs signals*/
		int    *data  = new int[m*nsymb];			//Coded bits bits
		double *mod   = new double[2 * (nsymb)];	//Constellation points
		double *modwp = new double[2 * (npsymb)];	//Constellation points with pilots
		int* indexes  = new int[npsymb];

		double *tx = new double[2 * (int)nofsamp]; // Transmitted samples
		double *rx = new double[2 * (int)nofsamp]; //RX samples (double buffer)

		double *modwpRX1 = new double[2 * npsymb];	// Output of timing
		double *modwpRX2 = new double[2 * npsymb];	// Output of Frame sync and Frequency offset
		double *modwpRX  = new double[2 * npsymb];	// Output of equalizer
		double *modRX    = new double[2 * nsymb];	// Output of pilot stripping
		int    *dataRX   = new int[m*nsymb];		// Estimated coded bits


		bool cont = true;
		int specmeas = 0;
		int f,i;
		double pow = 1.;
		int tt = 0;
		double ph1 = 0.;
		int fp=0;
		int nf;
		int dxt=0, dxr=0;

		pt1 = 0;
		pt2 = sfifo / 2;
		pr1 = 0;
		pr2 = sfifo / 2;

		std::thread *t1;
		if (actTX == 1 && actRX == 1)
		{
			delete[] TX; TX = (double*)calloc(2 * nofsamp*sfifo,8);
			delete[] RX; RX = (double*)calloc(2 * nofsamp*sfifo, 8);
			USRP->Initiate();		
			t1=new thread(TXRX);
		}
		else if (actTX == 1 && actRX == 0)
		{
			delete[] TX; TX = (double*)calloc(2 * nofsamp*sfifo,8);
			t1=new thread(Send);
		}
		else if (actRX == 1&& actTX == 0)
		{
			delete[] RX; RX = (double*)calloc(2 * nofsamp*sfifo, 8);
			USRP->Initiate();
			t1=new thread(Receive);
		}
		else
			continue;

		__int64 avewaittx = 0;
		__int64 avewaitrx = 0;
		int losttx = 0;
		int lostrx = 0;
		int frec;
		int ndec = 0;
		double RXpower = 1;
		int nerr = 0;
		FILE* ft;
		for(f=1;f<fmax && cont;f++)
		{
			speed.start();

			// Keyboard interaction
			if (f % 10 == 0)// Interaction with keyboard
			{
				if (nerr != USRP->nerr)
				{
					printf("USB: Lost packet (%d,%d)\n",USRP->nerr,f);
					nerr = USRP->nerr;
				}
				while (_kbhit())
				{
					switch (_getch())
					{
					case '+':if (actTX) {USRP->SetGain(1, 0, ++gainT); gainT = USRP->GetGain(1, 0);printf("GainT=\t%f\n", gainT); } break;
					case '-':if (actTX) {USRP->SetGain(1, 0, --gainT); gainT = USRP->GetGain(1, 0);printf("GainT=\t%f\n", gainT); } break;
					case ']':if (actRX) {USRP->SetGain(0, 0, ++gainR); gainR = USRP->GetGain(0, 0);printf("GainR=\t%f\n", gainR); }break;
					case '[':if (actRX) {USRP->SetGain(0, 0, --gainR); gainR = USRP->GetGain(0, 0);printf("GainR=\t%f\n", gainR); }break;


					case 'b':
						if (actRX == 0)break;
						BER[0].Display();
						break;

					case 'B': // Reset BER meter
						if (actRX == 0)break;
						printf("Reset BER\n");
						BER[0].Reset();
						break;

					case 'P': // Store current BER statistics
						BER[0].Display();
						printf("Append BER stat in BER_results.txt\n");
						ft = fopen("Data/BER_stat.txt", "a");
						a.Print_Parameters_on_Line(ft);
						BER[0].Display_on_Line(ft);
						fprintf(ft, "\n");
						fclose(ft);
						break;
					case 'R': // Reset receiver
						if (actRX == 0)break;
						printf("Reset Receiver\n");
						FSYNC_FOFF->Reset();
						FSYNC_FOFF->pointer += Nd / 2;
						CT->Reset();
						MMSE_PLL->Reset();
						break;					
					case 'm':
						specmeas ^= 1;
						if (specmeas == 1)printf("Now measuring spectra.\n");
						else
						{
							printf("Now NOT measuring spectra.\n");
							if (actTX){SP[0].Print(); SP[0].Reset();}
							if (actRX)
							{
								SP[1].Print();SP[1].Reset();
								SP[2].Print();SP[2].Reset();					
							}
						}
						break;

					case 'u': // Output
						fprintf(stdout,"USRP: Total errors: \t%d\n", USRP->nerr);
						if (actTX)
						{
							fprintf(stdout, "USRP: TX errors   : \t%d\n", USRP->nerrTXOF);
							printf("PC   TX:\t%d\nUSRP TX:\t%d\t(D=%d)\n", pt1, pt2, dxt);
							printf("Wait average at TX=\t%f\n", (double)avewaittx / pt1);
							printf("Lost at TX (PC)   =\t%d\n", losttx);
							printf("GainT             =\t%f\n", gainT);
						}
						if (actRX)
						{
							fprintf(stdout, "USRP: RX errors   : \t%d\n", USRP->nerrRXOF);
							printf("PC   RX:\t%d\nUSRP RX:\t%d\t(D=%d)\n", pr1, pr2, dxr);
							printf("Wait average at RX=\t%f\n",(double) avewaitrx/pr1);
							printf("Lost at RX (PC)   =\t%d\n", lostrx);
							printf("Received frame    =\t%d\n", frec);
							printf("GainR             =\t%f\n", gainR);
						}
						break;

					case 'p': // Output
						printf("Store waveforms in data files.\n");
						fp = 1;
						break;
					case 'r': // Rewind output
						if (actRX == 0)break;
						fprintf(stdout, "Reinitialize IQ plot and Scattering diagram files.\n");
						fclose(out1); fclose(out2); fclose(out3); fclose(out4);
						out1 = fopen("Data/IQplot.txt", "w");
						out2 = fopen("Data/scatdg_BEFORE.txt", "w");
						out3 = fopen("Data/scatdg_MMSE.txt", "w");
						out4 = fopen("Data/Centroids.txt", "w");
						break;

					case 'q':
						printf("Quit current scenario.\n");
						cont = false;
						break;
					case 'U':
						if (actRX == 0)break;
						printf("Update centroid locations into demodulator.\n");
						CT->Set_Demodulator(RXMOD, DEM);
						break;
					case 'S': // Output
						a.Print_Parameters(); break;

					case 's': // Output
						if (actRX == 0)break;
						printf("%d: Status of receiver______________________________________\n", f);
						AGC			->Display();
						TSYNC		->Display();
						FSYNC_FOFF	->Display();
						MMSE_PLL	->Display();
						CT			->Display_Centroids(stdout, true);
						break;
					case ' ':
						speed.Print(nofsamp/1000., "ksamp/sec");
						printf("_______ Available keys______________________\n");
						printf("-+: Decrease/Increase USRP TX gain.\n");
						printf("[]: Decrease/Increase USRP RX gain.\n");
						printf("S : Displays input parameters.\n");
						printf("s : Displays receiver status.\n");
						printf("R : Reset receiver\n");
						printf("b : Displays BER statistics.\n");
						printf("B : Reset BER statistics.\n");
						printf("P : Append BER stat in BER_results.txt.\n");
						printf("p : Stores IQ plot and Scattering diagram.\n");
						printf("r : Reinitializes IQ and scatdg files.\n");
						printf("U : Update centroid locations into demodulator.\n");
						printf("u : Status of USRP and PC TX and RX.\n");
						printf("m : Toggle spectrum mesurements.\n");
						printf("q : Quit current scenario.\n");
						printf("____________________________________________\n");
						break;
					default: break;
					}
				}
			}

			//Transmitter (TX)
			if (actTX > 0)
			{
				for (i = 0; i < 32; i++)data[i] = (f >> i) & 1;  // Encode frame number
				Source  ->seed = f;								// Set Seed to frame number
				Source	->Run(m*nsymb-32, data+32);					// Source
				MOD		->Run(nsymb, data, mod);				// Modulation
				PILOT   ->Run(-npsymb, mod, modwp);				// Pilot insertion


				// Sync with USRP
				SyncTX(avewaittx, losttx);
				TXfil	->RunImpulse(npsymb, modwp, TX+(pt1&(sfifo-1))*2*nofsamp, ns);	// Shaping filter
				if (specmeas)SP[0].Run(nofsamp, TX + (pt1&(sfifo - 1)) * 2 * nofsamp);
				pt1++;		
			}

			//Receiver (RX)
			if (actRX > 0)
			{
				if (actTX == 0)SyncRX(avewaitrx, lostrx);
				RXfil->Run(nofsamp, RX+(pr1%sfifo)*2*nofsamp, rx);	// SRRC filter
				pr1++;
				if (specmeas)
				{
					SP[1].Run(nofsamp, RX + (pr1%sfifo) * 2 * nofsamp);
					SP[2].Run(nofsamp, rx);
				}
	//			continue;
				AGC->Run(nofsamp, rx, rx);
				nsx = TSYNC->Run(nofsamp, rx, modwpRX1);			// Timing
				nf  = FSYNC_FOFF->Run(nsx, modwpRX1,modwpRX2);		// Frame and coarse foffset recovery
				if(FSYNC_FOFF->ini==1)continue;
				if (nf == 0)continue;								// No frames available
				MMSE_PLL->Run(npsymb, modwpRX2, modwpRX,indexes);   // MMSE and PLL, returns also the indexes of transmitted symbols for CT
				CT      ->Update_Positions(npsymb, modwpRX, indexes);     // Update Centroids
				if (fp > 0)
				{
					for (i = 0; i < 200; i++)
					{
						fprintf(out1, "%lf\t%lf\n", rx[2 * i], rx[2 * i + 1]);
						fprintf(out2, "%lf\t%lf\n", modwpRX2[2 * i], modwpRX2[2 * i + 1]);
						fprintf(out3, "%lf\t%lf\n", modwpRX[2 * i], modwpRX[2 * i + 1]);
					}
					if (fp == 1)
					{
						CT->Display_Centroids(out4, true);
					}
					fflush(out1);
					fflush(out2);
					fflush(out3);
					fflush(out4);
					fp--;
				}

				PILOT	->RunInv(-npsymb, modwpRX, modRX);			// Pilot stripping
				DEM		->RunSoft(nsymb, modRX, dataRX, true);			// Soft Demodulator
				for (i = frec =0; i < 32; i++)
				{
					if (dataRX[i] > 0)
					{
						frec += 1 << i;
						data[i] = 1;
					}
					else
					{
						data[i] = 0;
					}
				}
				Source->seed = frec;
				Source->Run(m*nsymb - 32, data + 32);					// Source
				BER[0].Run(m*nsymb,data,dataRX);				// BER Meter
			}
		}
		mtx.lock();
		pt2 = pr2 = -100;
		mtx.unlock();
		a.Print_Parameters_on_Line(file);
		if (actRX > 0)
		{
			printf("%d: Status of receiver______________________________________\n", f);
			TSYNC->Display();
			FSYNC_FOFF->Display();
			MMSE_PLL->Display();
			CT->Display_Centroids(stdout, true);
			BER[0].Display();
			fprintf(stdout, "USRP: RX errors   : \t%d\n", USRP->nerrRXOF);
			printf("PC   RX:\t%d\nUSRP RX:\t%d\t(D=%d)\n", pr1, pr2, dxr);
			printf("Wait average at RX=\t%f\n", (double)avewaitrx / pr1);
			printf("Lost at RX        =\t%d\n", lostrx);
			printf("Received frame    =\t%d\n", frec);
			printf("GainR             =\t%f\n", gainR);
		}
		if (actTX)
		{
			fprintf(stdout, "USRP: TX errors   : \t%d\n", USRP->nerrTXOF);
			fprintf(stdout, "USRP: TX errors   : \t%d\n", USRP->nerrTXOF);
			printf("PC   TX:\t%d\nUSRP TX:\t%d\t(D=%d)\n", pt1, pt2, dxt);
			printf("Wait average at TX=\t%f\n", (double)avewaittx / pt1);
			printf("Lost at TX        =\t%d\n", losttx);
			printf("GainT             =\t%f\n", gainT);
		}
		BER[0].Display_on_Line(file);
		fprintf(file,"%f\t", gainT);
		fprintf(file,"%f\t", gainR);
		fprintf(file, "%d\t", losttx);
		fprintf(file, "%d\t", lostrx);
		fprintf(file, "%d\t", USRP->nerrTXOF);
		fprintf(file, "%d\t", USRP->nerrRXOF);
		fprintf(file, "\n");
		fflush(file);


		delete[] data, mod, modwp, tx;
		delete[] dataRX, modRX, modwpRX, modwpRX1, modwpRX2, rx;
		delete[] indexes;
	
		delete Source, MOD, PILOT, TXfil;
		delete RXfil,TSYNC,FSYNC_FOFF,MMSE_PLL;


		delete[] SP;
		delete[] BER;
		USRP->Close();
		delete USRP;	
	}while(a.Update_Parameters());
	fclose(out3);
	fclose(out2);
	fclose(out1);
	fclose(file);
}