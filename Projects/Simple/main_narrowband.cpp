// Include TOPCOM++ headers of used blocks
#include "Parameter_Manager.h"
#include "PN_Source.h"
#include "Modulator.h"
#include "AWGN_Channel.h"
#include "Demodulator.h"
#include "BER_meter.h"
#include "Filter.h"
#include "Spectrum.h"
#include "speed.h"
int main()
{
	/* Declaration of global variables								   */
	/***********************************************************/
	double ebn0db,rolloff;
	int m,Nfil,ns,i;
	int specmeas = 0;
	char nameout[132];

	double nuf, nua, A;
	double tta, ttf;
	double abs;
	double ch[2];
	const double tpi = 6.283185307179586476925286766559;
	double temp;



	/* Compile-time specified parameters					   */
	/***********************************************************/
	int N=500;  // Frame size

	/* Instantiation of blocks NOT affected by users parameter */
	/***********************************************************/

	/* Allocation of buffers NOT affected by users parameter */
	/***********************************************************/

	/* Parameter specified by the user (using Parameter_Manager)*/
	/***********************************************************/
	Parameter_Manager a(40);
	FILE *input;
	char nameinput[132] = "Data/input_narrowband.txt";
	input = fopen(nameinput, "r+");
	if (input == NULL)input = fopen(nameinput, "w+");

	 a.Add_Parameter_From_File(input,"Name output file",nameout);	
	 a.Add_Parameter_From_File(input,"Eb/N0 [dB]",&ebn0db);
	 a.Add_Parameter_From_File(input,"Number of modulation bits", &m);
	 a.Add_Parameter_From_File(input,"SRRC roll-off", &rolloff);
	 a.Add_Parameter_From_File(input,"Length of filters", &Nfil);
	 a.Add_Parameter_From_File(input,"Number of samples per symbol", &ns);
	 a.Add_Parameter_From_File(input, "Normalized Frequency offset", &nuf);
	 a.Add_Parameter_From_File(input, "Amplitude fluctuations frequency (nua)", &nua);
	 a.Add_Parameter_From_File(input, "Amplitude PP fluctuations", &A);

	 // ...
	fclose(input);
	if (a.error)
	{
		printf("Main: The file %s was not ready. Fill up all the required parameters.\n", nameinput);
		return 1;
	}

	FILE* file = fopen(nameout, "w");

	/***********************************************************/
	/* Parameters loop  */
	/***********************************************************/
	do // 
	{
		FILE* out = fopen("Data/IQplot.txt", "w");
		FILE* out2 = fopen("Data/scatdg.txt", "w");
		a.Print_Parameters();
		if (ns <= 1)continue;

		/* Compute derived parameters */
		/***********************************************************/

		/* Instantiation of blocks affected by users parameter    */
		/***********************************************************/
		PN_Source* Source	= new PN_Source;

		Modulator* Mod	= QAM_Modulator(m);

		Filter* TXFil = 0;
		Filter* RXFil = 0;
		TXFil = SRRC(rolloff, ns, Nfil, true);
		TXFil->Set_Unitary_Energy();

		RXFil = SRRC(rolloff, ns, Nfil, true);
		RXFil->Set_Unitary_Energy();

		AWGN_Channel* AWGN	= new AWGN_Channel;
		AWGN->Set_EsN0dB(ebn0db+10.*log10((double)m));

		Demodulator* Demod	= new Demodulator;
		Demod->SetParameters(Mod);

		BER_meter* BER	= new BER_meter;
		int delay = 0;
		if(ns>1)delay = 2 * Nfil*m;
		BER->SetParameters(delay, 30);

		Spectrum* SP = new Spectrum[3];
		SP[0].SetParameters(12, 1. / ns, "Data/spectrumTX.txt");
		SP[1].SetParameters(12, 1. / ns, "Data/spectrumRX.txt");
		SP[2].SetParameters(12, 1. / ns, "Data/spectrumRX_afterSRRC.txt");
		specmeas = 0;


		Speed speed;

		/* Allocation of buffers affected by users parameter	   */
		/***********************************************************/
		int* data	= new int	 [N*m];
		double* mod = new double[2 * N];
		double* tx = new double[2 * N*ns];
		double* rx = new double[2 * N*ns];
		double* rx2 = new double[2 * N*ns];
		double* modRX = new double[2 * N];
		int* dataRX	= new int	 [N*m];

		/* Initializations										   */
		/***********************************************************/

		/***********************************************************/
		/* Simulation loop										   */
		/***********************************************************/
		int f;
		bool cont=true;
		tta = ttf = 0.;

		for (f = 0;cont;f++)
		{

			if (f % 10 == 0)// Interaction with keyboard
			{
				if (_kbhit())
				{
					char cc = _getch();
					switch (cc)
					{
					case 'b':BER->Display(); break;

					case 'p':
						fprintf(stdout, "Storing waveform of 200 symbols in IQplot.txt\n");
						fprintf(out, "\n\n");
						for (i = 0; i < 200 * ns; i++)
						{
							if (i % (2 * ns) == 0)
							{
								fprintf(out, "%f\t%f\t%f\n", (double)2.*ns/ns, rx2[2 * i], rx2[2 * i + 1]);
								fprintf(out, "\n");
								fprintf(out, "%f\t%f\t%f\n", 0., rx2[2 * i], rx2[2 * i + 1]);
							}
							else
							{
							    fprintf(out, "%f\t%f\t%f\n", (double)(i%(2*ns))/ns, rx2[2 * i], rx2[2 * i + 1]);
							}
						}
						fflush(out);
						fprintf(stdout, "Storing scattering diagram of 200 symbols in scatdg.txt\n");
						fprintf(out2, "\n\n");
						for (i = 0; i < 200; i++)
						{
							fprintf(out2, "%f\t%f\n", modRX[2 * i], modRX[2 * i + 1]);
						}
						fflush(out2);
						break;
					case 'r': // Rewind output
						fprintf(stdout, "Reinitialize data files.\n");
						fclose(out); fclose(out2); 
						out = fopen("Data/IQplot.txt", "w");
						out2 = fopen("Data/scatdg.txt", "w");
						break;

					case 'n':
						cont = false;
						break;

					case 'm':
						specmeas ^= 1;
						if (specmeas)fprintf(stdout, "Spectrum measurements turned ON.\n");
						else        fprintf(stdout,  "Spectrum measurements turned OFF.\n");
						break;
	
					default:
						speed.Print((double)N*m / 1000, "kbit/sec");
						fprintf(stdout, "______________Available keys__________________\n");
						fprintf(stdout, "'r' Reset IQ plot and and scattering diagram. .\n");
						fprintf(stdout, "'b' Displays BER statistic.\n");
						fprintf(stdout, "'p' store IQ plot samples and scattering diagram.\n");
						fprintf(stdout, "'m' Toggle spectrum measures.\n");
						fprintf(stdout, "'n' Skip to next scenario.\n");
						break;
					}
					_getch();
				}
			}

			// TX
			Source	->Run(N*m, data);
			Mod->Run(N, data, mod);
			TXFil->RunImpulse(N, mod, tx,ns);
			if(specmeas)SP[0].Run(N*ns, tx);

			// Narrowband sinusoidal channel model 
			for (i = 0; i < N*ns; i++)
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

			// Channel
			AWGN	->Run(2*N*ns, tx, rx);

			// RX
			if(specmeas)SP[1].Run(N*ns, rx);

			RXFil->Run(N*ns, rx, rx2);

			if(specmeas)SP[2].Run(N*ns, rx2);

			// Decimation
			for (i = 0; i < N; i++)
			{
				modRX[2 * i] = rx2[2 * (i*ns)];
				modRX[2 * i + 1] = rx2[2 * (i*ns) + 1];
			}

			Demod	->Run(N, modRX, dataRX);
			BER->Run(N*m, data, dataRX);
//			if (BER->IsReliable())break;
		}

		/* Post processing single run							   */
		/***********************************************************/
		BER->Display();

		a.Print_Parameters_on_Line(file);
		BER->Display_on_Line(file);
		fprintf(file, "\n");
		fflush(file);
		fclose(out);
		fclose(out2);


		/* Delete blocks and buffers affected by users parameter   */
		/***********************************************************/
		delete Source;
		delete Mod;
		delete AWGN;
		delete TXFil;
		delete RXFil;
		delete Demod;
		delete BER;

		delete[] SP;
		delete[] data;
		delete[] mod;
		delete[] tx;

		delete[] rx;
		delete[] rx2;
		delete[] modRX;
		delete[] dataRX;
		//...
	} while (a.Update_Parameters());

	/* Global Post processing								   */
	/***********************************************************/

	/* Delete blocks and buffers NOT affected by users parameter   */
	/***********************************************************/
	return 0;
}
