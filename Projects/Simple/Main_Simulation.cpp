#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
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
#include "FIFO.h"
#include "Pilot.h"
#include "Equalizer.h"
#include "Frame_sync.h"
#include "Automatic_Gain_Control.h"

int main()
{
	int N, D, H, P;
	int fmax;
	double nua, nuf, A, tta, ttf, ch[2], abs, temp;
	const double tpi = 6.283185307179586476925286766559;
	char nameinput[128] = "Data/input_Simulation.txt";
	double alpha;
	double decw = 0.; // If different it alters the power computation
	double alphap = 0.1;// 0.005;
	double decwp = 0.5;//0.8;
	double afoff = 0.05;
	int coh = 5;
	bool specmeas = false;// Turn this on to turn on spectrum measures
	int Nd = 1;
	int Nfil = 5;
	int L0;
	double rolloff;
	double bwav = 1e-2;// Bandwidth for performins averages
	int m, ns, nsymb, npsymb, nsx;
	int   nofsamp;

	double SNR,ebn0;
	int ncand = 5;
	char name[132];
	double gamma;

	Parameter_Manager a(40);
	FILE *input;
	input = fopen(nameinput, "r+");
	if (input == NULL)input = fopen(nameinput, "w+");
	a.Add_Parameter_From_File(input, "Name of output file", name);
	a.Add_Parameter_From_File(input, "Upsampling factor", &ns);
	a.Add_Parameter_From_File(input, "Number of frames", &fmax);
	a.Add_Parameter_From_File(input, "Number of data sections in a frame (N)", &N);
	a.Add_Parameter_From_File(input, "Header length (H)", &H);
	a.Add_Parameter_From_File(input, "Data Length (D)", &D);
	a.Add_Parameter_From_File(input, "Pilot length (P)", &P);
	a.Add_Parameter_From_File(input, "Eb/N0 [dB]", &ebn0);

	a.Add_Parameter_From_File(input, "Normalized Frequency offset", &nuf);
	a.Add_Parameter_From_File(input, "Amplitude fluctuations normalized frequency", &nua);
	a.Add_Parameter_From_File(input, "Amplitude PP fluctuations", &A);
	a.Add_Parameter_From_File(input, "SRRC filter roll-off", &rolloff);
	a.Add_Parameter_From_File(input, "QAM modulation efficiency", &m);
	a.Add_Parameter_From_File(input, "Upsampling factor", &ns);
	a.Add_Parameter_From_File(input, "AGC bandwidth", &gamma);
	a.Add_Parameter_From_File(input, "Observation window of timing", &L0);
	a.Add_Parameter_From_File(input, "Updating step of frequency estimator", &afoff);
	a.Add_Parameter_From_File(input, "Number of coherent sums (nc)", &coh);
	a.Add_Parameter_From_File(input, "Updating step of Equalizer", &alpha);
	a.Add_Parameter_From_File(input, "Updating step of PLL", &alphap);
	a.Add_Parameter_From_File(input, "Weight of decision aided error for Equalizer", &decw);
	a.Add_Parameter_From_File(input, "Weight of decision aided error for", &decwp);
	fclose(input);
	if (a.error)
	{
		printf("Main: The file %s was not ready. Fill up all the required parameters.\n", nameinput);
		exit(1);
	}
	FILE* file = fopen(name, "w");
	FILE* out1 = fopen("Data/IQplot.txt", "w");
	FILE* out2 = fopen("Data/scatdg_BEFORE.txt", "w");
	FILE* out3 = fopen("Data/scatdg_MMSE.txt", "w");
	a.Print_Parameters_on_Line(file, true);
	fprintf(file, "\n");
	do {


		if (Nd % 2 != 1)continue;
		if (coh < 1)continue;
		if (coh > H / 2)continue;
		if (nuf == 0.)coh = H / 2;
		else        coh = (int)(1 / fabs(nuf) / 10);
		if (coh > H / 2)coh = H / 2;

		a.Print_Parameters();
		nsymb = N * D;
		npsymb = H + N * (D + P);
		nofsamp = ns * npsymb;

		Speed speed;

		PN_Source *Source;

		// Modulation
		Modulator *MOD = QAM_Modulator(m);

		// Pilot Insertion
		Pilot* PILOT = new Pilot;
		PILOT->SetRegularFrame(N, H, D, P); // Frame structure with header of length H and N repetitions of D+P
		int    *dd = new int[m*(H + N * P)]; //source bits
		double* pilot = new double[2 * (H + N * P)];
		Source = new PN_Source((int)(log((H + N * P)*m) / log(2.) + 0.9999));
		Source->Run((H + N * P)*m, dd);
		MOD->Run((H + N * P), dd, pilot);
		PILOT->SetPilot((H + N * P), pilot);
		delete[] dd;
		delete[] pilot;
		delete Source;

		Source = new PN_Source;

		Filter* TXfil = SRRC(rolloff, ns, Nfil, true);
		TXfil->Set_Unitary_Energy();
		SNR = ebn0 + 10 * log10((double)m);
		AWGN_Channel* AWGN = 0;
		AWGN = new AWGN_Channel;
		AWGN->Set_EsN0dB(SNR);

		Filter* RXfil = SRRC(rolloff, ns, Nfil, true);
		RXfil->Set_Unitary_Energy();
		RXfil->SetMatched();

		// AGC
		Automatic_Gain_Control* AGC = new Automatic_Gain_Control;
		AGC->SetParameters(1, gamma);

		// Timing
		Time_Synchronizer* TSYNC = new Time_Synchronizer;
		TSYNC->SetParameters(ns, L0);

		// Frame Sync and frequency offset estimation
		Frame_Sync* FSYNC_FOFF = new Frame_Sync;
		double* ppp = PILOT->GetPilotSequence(H + N * (D + P));
		FSYNC_FOFF->SetParameters(npsymb, H, ppp, ncand, coh, afoff);

		// MMSE_PLL 
		Equalizer* MMSE_PLL = new Equalizer;
		MMSE_PLL->SetParameters(MOD, Nd, alpha, 1, 0, 0., decw);
		MMSE_PLL->SetTraining(npsymb, ppp);
		MMSE_PLL->ActivatePLL(alphap, decwp);
		delete[] ppp;
	
		// Demodulator
		Demodulator* DEM = new Demodulator;
		Modulator* RXMOD = new Modulator;
		RXMOD->SetParameters(MOD->m, MOD->constellation, MOD->mapping);
		DEM->SetParameters(RXMOD);
		DEM->SetSigma(pow(10., -SNR / 10) / 2.);

		Delay* DEL = new Delay[1];
		DEL[0].SetParameters(N*m*D, 0);

		BER_meter* BER = new BER_meter[1];
		BER[0].SetParameters(0, 30);
		BER[0].SetFrameSize(N*D*m);
		BER[0].SetSoft();


		Spectrum* SP = new Spectrum[3];
		SP[0].SetParameters(10, 1. / ns, "Data/spectrumTX.txt");
		SP[1].SetParameters(10, 1. / ns, "Data/spectrumRX.txt");
		SP[2].SetParameters(10, 1. / ns, "Data/spectrumRX_afterSRRC.txt");
		specmeas = 0;


		/* Allocate buffers for storing input and outputs signals*/
		int    *data = new int[m*nsymb]; //source bits
		double *mod = new double[2 * (nsymb)];	 //Constellation points
		double *modwp = new double[2 * (npsymb)];	//Constellation points

		double *tx = new double[2 * (int)nofsamp]; // Transmitted samples
		double *rx = new double[2 * (int)nofsamp]; //RX samples (double buffer)

		double *modwpRX2 = new double[2 * npsymb];	// Output of timing
		double *modwpRX1 = new double[2 * npsymb];	// Output of timing
		double *modwpRX = new double[2 * npsymb];	// Output of equalizer
		double *modRX = new double[2 * nsymb];	// Output of pilot stripping
		int    *dataRX = new int[m*nsymb];			// Estimated coded bits


		bool cont = true;
		int f, i;
		double pow = 1.;
		int tt = 0;
		double ph1 = 0.;
		int fp = 0;
		int nf;
		tta = ttf = 0.;
		bool decoding = false;

		double ma = 1., mt = 0., mf = 0., mth = 0.;
		double sa = 1., st = 0., sf = 0., sth = 0.;
		int frec;
		int fx = 0;
		for (f = 1; f <= fmax && cont; f++)
		{
			speed.start();

			// Keyboard interaction
			if (f % 10 == 0)// Interaction with keyboard
			{
				while (_kbhit())
				{
					switch (_getch())
					{
					case 'R': // Reset BER meter
						printf("Reset BER statistics\n");
						BER->Reset();
						break;
					case 'b':
						BER[0].Display();
						break;
					case 'p': // Output
						printf("Store waveforms in data files.\n");
						fp = 1;
						break;
					case 'd':
						decoding ^= 1;
						if (decoding == 1)printf("Now decoding \n");
						else              printf("Now NOT decoding \n");
						break;

					case 'm':
						specmeas ^= 1;
						if (specmeas)fprintf(stdout, "Spectrum measurements turned ON.\n");
						else        fprintf(stdout, "Spectrum measurements turned OFF.\n");
						break;

					case 'r': // Rewind output
						fprintf(stdout, "Reinitialize data files.\n");
						fclose(out1); fclose(out2); fclose(out3); 
						out1 = fopen("Data/IQplot.txt", "w");
						out2 = fopen("Data/scatdg_BEFORE.txt", "w");
						out3 = fopen("Data/scatdg_MMSE.txt", "w");
						break;

					case 'q':
						printf("Quit current scenario.\n");
						cont = false;
						break;
					case 'S': // Output
						a.Print_Parameters();
						break;
					case 's': // Output
						printf("%d: Status of receiver______________________________________\n", f);
						AGC->Display();
						TSYNC->Display();
						FSYNC_FOFF->Display();
						MMSE_PLL->Display();
						break;
					case 'e': // Print estimates
						fprintf(stdout, "AGC  :\t%f\t%f\t%e\n", AGC->gain, ma, sqrt(sa - ma * ma));
						fprintf(stdout, "Tau  :\t%f\t%f\t%e\n", TSYNC->tau, mt, sqrt(st - mt * mt));
						fprintf(stdout, "Foff :\t%f\t%f\t%e\n", FSYNC_FOFF->foff / tpi, mf / tpi, sqrt(sf - mf * mf) / tpi);
						fprintf(stdout, "Theta:\t%f\t%f\t%e\n", MMSE_PLL->theta / tpi, mth / tpi, sqrt(sth - mth * mth) / tpi);
						break;
					default:
						speed.Print(npsymb*ns / 1000., "ksamp/sec");
						printf("_______ Available keys______________________\n");
						printf("S: Displays input parameters.\n");
						printf("s: Displays receiver status.\n");
						printf("b: Displays BER statistics.\n");
						printf("e: Displays Estimates and variances\n");
						printf("p: Stores waveforms in data files.\n");
						printf("m: Toggle spectrum measures.\n");
						printf("R: Reset BER statistics\n");
						printf("r: Reinitializes data files.\n");
						printf("q: Quit current scenario.\n");
						printf("____________________________________________\n");
						break;
					}
					_getch();
				}
			}

			//Transmitter (TX)
			Source->Run(m*nsymb, data);					// Source
			MOD->Run(nsymb, data, mod);				// Modulation
			PILOT->Run(-npsymb, mod, modwp);				// Pilot insertion
			TXfil->RunImpulse(npsymb, modwp, tx, ns);	// Shaping filter	
			if (specmeas)SP[0].Run(npsymb*ns, tx);

			DEL[0].Run(N*D*m, data, data);						// Delay Data	

			// Narrowband channel model 
			for (i = 0; i < npsymb*ns; i++)
			{
				abs = 1. + A * sin(tta*tpi);
				ch[0] = abs * cos(ttf*tpi);
				ch[1] = abs * sin(ttf*tpi);
				temp = ch[0] * tx[2 * i] - ch[1] * tx[2 * i + 1];
				tx[2 * i + 1] = ch[0] * tx[2 * i + 1] + ch[1] * tx[2 * i];
				tx[2 * i] = temp;
				tta = fmod(tta + nua / ns, 1);
				ttf = fmod(ttf + nuf / ns, 1.);
			}
			AWGN->Run(2 * npsymb*ns, tx, rx);

			if (specmeas)SP[1].Run(npsymb*ns, rx);

			//Receiver (RX)
			AGC->Run(npsymb*ns, rx, rx);						// Automatic Gain Control
			RXfil->Run(npsymb*ns, rx, rx);						// SRRC filter
			if (specmeas)SP[2].Run(npsymb*ns, rx);
			nsx = TSYNC->Run(npsymb*ns, rx, modwpRX1);			// Timing
			if (nsx != npsymb)printf("nsx=%d\n", nsx);
			int temp = FSYNC_FOFF->ini;
			nf = FSYNC_FOFF->Run(nsx, modwpRX1, modwpRX2);		// Frame and frequency offset recovery
			if (FSYNC_FOFF->ini == 1)continue;
			FSYNC_FOFF->Display();
			//if (temp^FSYNC_FOFF->ini)
			//{
			//	printf("Frame aquired in %d\n", f-fx);
			//	fx = f;
			//	FSYNC_FOFF->ini = 0;
			//	FSYNC_FOFF->Reset();
			//	continue;
			//}
			continue;
			if (nf == 0)continue;
			MMSE_PLL->Run(npsymb, modwpRX2, modwpRX);	// MMSE and PLL 
			PILOT->RunInv(-npsymb, modwpRX, modRX);		// Pilot stripping
			DEM->RunSoft(nsymb, modRX, dataRX, true);	// Soft Demodulator

			if (f < 4)continue;
			BER[0].Run(N*D*m, data, dataRX);					// BER Meter
			ma += bwav * (AGC->gain - ma);
			sa += bwav * (AGC->gain*AGC->gain - sa);

			mf += bwav * (FSYNC_FOFF->foff - mf);
			sf += bwav * (FSYNC_FOFF->foff*FSYNC_FOFF->foff - sf);

			mt += bwav * (TSYNC->tau - mt);
			st += bwav * (TSYNC->tau*TSYNC->tau - st);

			mth += bwav * (MMSE_PLL->theta - mth);
			sth += bwav * (MMSE_PLL->theta*MMSE_PLL->theta - sth);

			// Print IQ plot and scatdg	
			if (fp > 0)
			{
				for (i = 0; i < 500; i++)
				{
					fprintf(out1, "%lf\t%lf\n", rx[2 * i], rx[2 * i + 1]);
					fprintf(out2, "%f\t%f\n", modwpRX2[2 * i], modwpRX2[2 * i + 1]);
					fprintf(out3, "%f\t%f\n", modwpRX[2 * i], modwpRX[2 * i + 1]);
				}
				fflush(out1);
				fflush(out2);
				fflush(out3);
				fp--;

			}

			if (BER[0].IsReliable())break;
		}
		printf("%d: Status of receiver______________________________________\n", f);
		AGC->Display();
		TSYNC->Display();
		FSYNC_FOFF->Display();
		MMSE_PLL->Display();
		BER[0].Display();

		a.Print_Parameters_on_Line(file);
		fprintf(file, "%f\t%f\t%e\t", AGC->gain, ma, sqrt(sa - ma * ma));
		fprintf(file, "%f\t%f\t%e\t", TSYNC->tau, mt, sqrt(st - mt * mt));
		fprintf(file, "%f\t%f\t%e\t", FSYNC_FOFF->foff / tpi, mf / tpi, sqrt(sf - mf * mf) / tpi);
		fprintf(file, "%f\t%f\t%e\t", MMSE_PLL->theta / tpi, mth / tpi, sqrt(sth - mth * mth) / tpi);
		fprintf(file, "%f\t%f\t%e\t", MMSE_PLL->MSE);
		BER[0].Display_on_Line(file);
		fprintf(file, "\n");
		fflush(file);

		delete[] data, mod, modwp, tx;
		delete[] dataRX, modRX, modwpRX, modwpRX1, modwpRX2, rx;

		delete AWGN;

		delete Source,  MOD, PILOT, TXfil;
		delete RXfil, TSYNC, FSYNC_FOFF, MMSE_PLL;


		delete[] SP;
		delete[] BER;
		delete[] DEL;
	} while (a.Update_Parameters());
	fclose(out3);
	fclose(out2);
	fclose(out1);
	fclose(file);
}