#pragma once
#include "../NI-USRP/includes/niUSRP.h"
#include <stdio.h>
using namespace std;


class USRP_Manager
{

public:
	USRP_Manager(const int NT, const int NR, const double IQrate, const double carrierFrequency, const int nofsamp = 10000, const double gainT = 0, const double gainR = 30.);
	int TXRX(const short * dataBufferT, short * dataBufferR);
	int TX(const double * dataBufferT);
	int RX(double * dataBufferR);
	int RX(short * dataBufferR);
	int TXRX(const double* dataBufferT, double* dataBufferR);
	int SetGain(const int tx, const int ch, const double gain);
	double GetGain(const int tx, const int ch) const;
	int SetCarrier(const int tx, const double freq);
	int SetIQrate(const int tx, const int ch, const double IQrate);
	int Display(FILE * file=stdout);
	int Initiate();
	void Print_Attributes()
	{
		if (sessionHandleT != 0)
		{
			printf("TX 0 session:------------------------\n");
			Print_Attributes(sessionHandleT, "0");
			printf("TX 1 session:------------------------\n");
			Print_Attributes(sessionHandleT, "1");

		}
		if (sessionHandleR != 0)
		{
			printf("RX 0 session:------------------------\n");
			Print_Attributes(sessionHandleR, "0");
			printf("RX 1 session:------------------------\n");
			Print_Attributes(sessionHandleR, "1");
		}
	}


	int Close();
	int GetState(niUSRP_Session sessionHandle,niUSRP_ConstString channelList) const // Check if the driver is running
	{
		int temp;
		niUSRP_GetAttributeInt32(sessionHandle,channelList ,niUSRP_Attr_DriverState,&temp);
		return temp;
	}


	int nerr;
	int nerrTXOF; // Overflow internal buffer
	int nerrRXOF; // Overflow internal buffer
	int32_t numberOfDevicesInSystem;

private:
	int devtype;
	niUSRP_Session sessionHandleT;
	niUSRP_Session sessionHandleR;

	void Print_Attributes(niUSRP_Session sessionHandle, niUSRP_ConstString channelList);

	niUSRP_Device* devices;  //NI device

	void ManageError(niUSRP_Session sessionHandle);
	double coercedIQT0;
	double coercedIqRateT0;
	double coercedCarrierFrequencyT0;
	double coercedGainT0;

	double coercedIQT1;
	double coercedIqRateT1;
	double coercedCarrierFrequencyT1;
	double coercedGainT1;

	double coercedIQR0;
	double coercedIqRateR0;
	double coercedCarrierFrequencyR0;
	double coercedGainR0;

	double coercedIQR1;
	double coercedIqRateR1;
	double coercedCarrierFrequencyR1;
	double coercedGainR1;

	int nofsamp;
	int timeout;
	int NT,NR;
	double* buffer; // Store values
	niUSRP_String errorDescriptionBuffer;

};

