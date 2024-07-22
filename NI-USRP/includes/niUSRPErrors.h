/*
 * Errors for the NI-USRP C API
 *
 * Copyright (c) 2011-
 * National Instruments Corporation.
 * All rights reserved.
 */

#ifndef ___niUSRPErrors_h___
#define ___niUSRPErrors_h___


/****************************************************************************
 *  Error Codes 
 ****************************************************************************/

#define niUSRP_Error_Base                             0xBFFA4000
#define niUSRP_Status_Success                         (0)
                                                      
#define niUSRP_Error_InternalSoftwareFault            (niUSRP_Error_Base+1)
#define niUSRP_Error_MemoryAllocation                 (niUSRP_Error_Base+2)
#define niUSRP_Error_DriverLoadFailure                (niUSRP_Error_Base+3)
#define niUSRP_Error_NullParameter                    (niUSRP_Error_Base+4)
#define niUSRP_Error_FunctionNotSupported             (niUSRP_Error_Base+5)
#define niUSRP_Error_TimeoutExceeded                  (niUSRP_Error_Base+6)
#define niUSRP_Error_LateCommand                      (niUSRP_Error_Base+7)
#define niUSRP_Error_BrokenChain                      (niUSRP_Error_Base+8)
#define niUSRP_Error_Overflow                         (niUSRP_Error_Base+9)
#define niUSRP_Error_BadPacket                        (niUSRP_Error_Base+10)
#define niUSRP_Error_Underflow                        (niUSRP_Error_Base+11)
#define niUSRP_Error_SeqError                         (niUSRP_Error_Base+12)
#define niUSRP_Error_TimeError                        (niUSRP_Error_Base+13)
#define niUSRP_Error_PacketUnderflow                  (niUSRP_Error_Base+14)
#define niUSRP_Error_BurstSeqError                    (niUSRP_Error_Base+15)
#define niUSRP_Error_InvalidSessionHandle             (niUSRP_Error_Base+16)
#define niUSRP_Error_InvalidAttributeID               (niUSRP_Error_Base+17)
#define niUSRP_Error_InvalidAttributeScope            (niUSRP_Error_Base+18)
#define niUSRP_Error_UnsupportedDevice                (niUSRP_Error_Base+19)
#define niUSRP_Error_InvalidTimeSpecifier             (niUSRP_Error_Base+20)
#define niUSRP_Error_InvalidTriggerType               (niUSRP_Error_Base+21)
#define niUSRP_Error_InvalidSourceTerminal            (niUSRP_Error_Base+22)
#define niUSRP_Error_InvalidPolarity                  (niUSRP_Error_Base+23)
#define niUSRP_Error_DriverInRunningState             (niUSRP_Error_Base+24)
#define niUSRP_Error_DriverNotInRunningState          (niUSRP_Error_Base+25)
#define niUSRP_Error_MultipleAttributeGetConflict     (niUSRP_Error_Base+26)
#define niUSRP_Error_InvalidChannelName               (niUSRP_Error_Base+27)
#define niUSRP_Error_InvalidFetchSize                 (niUSRP_Error_Base+28)
#define niUSRP_Error_InternalDriverError              (niUSRP_Error_Base+29)
#define niUSRP_Error_LocalOscillatorFailedToLock      (niUSRP_Error_Base+30)
#define niUSRP_Error_ReferenceClockFailedToLock       (niUSRP_Error_Base+31)
#define niUSRP_Error_MimoLinkFailedToLock             (niUSRP_Error_Base+32)
#define niUSRP_Error_InvalidTimeout                   (niUSRP_Error_Base+33)
#define niUSRP_Error_InvalidSampleWidth               (niUSRP_Error_Base+34)
#define niUSRP_Error_InvalidHostDataType              (niUSRP_Error_Base+35)
#define niUSRP_Error_HostDataTypeConflict             (niUSRP_Error_Base+36)
#define niUSRP_Error_InvalidSizeOfDataBuffer          (niUSRP_Error_Base+37)
#define niUSRP_Error_Alignment                        (niUSRP_Error_Base+38)
#define niUSRP_Error_PpsWithGPSDC                     (niUSRP_Error_Base+39)
#define niUSRP_Error_GPSFailedToLock                  (niUSRP_Error_Base+40)
#define niUSRP_Error_UnsupportedLockSensor            (niUSRP_Error_Base+41)
#define niUSRP_Error_InvalidDeviceName                (niUSRP_Error_Base+42)
#define niUSRP_Error_InvalidOptionString              (niUSRP_Error_Base+43)
#define niUSRP_Error_InvalidSignal                    (niUSRP_Error_Base+44)
#define niUSRP_Error_IQRateDiffersAmongChannels       (niUSRP_Error_Base+45)
#define niUSRP_Error_OutOfSequence                    (niUSRP_Error_Base+46)
#define niUSRP_Error_UserPayload                      (niUSRP_Error_Base+47)
#define niUSRP_Error_Invalid290xDeviceName            (niUSRP_Error_Base+48)

#define niUSRP_Warning_InternalSoftwareFault          (-(int32_t)niUSRP_Error_InternalSoftwareFault)
#define niUSRP_Warning_TimeoutExceeded                (-(int32_t)niUSRP_Error_TimeoutExceeded)
#define niUSRP_Warning_Overflow                       (-(int32_t)niUSRP_Error_Overflow)
#define niUSRP_Warning_Underflow                      (-(int32_t)niUSRP_Error_Underflow)
#define niUSRP_Warning_SeqError                       (-(int32_t)niUSRP_Error_SeqError)
#define niUSRP_Warning_TimeError                      (-(int32_t)niUSRP_Error_TimeError)
#define niUSRP_Warning_PacketUnderflow                (-(int32_t)niUSRP_Error_PacketUnderflow)
#define niUSRP_Warning_BurstSeqError                  (-(int32_t)niUSRP_Error_BurstSeqError)
#define niUSRP_Warning_UserPayload                    (-(int32_t)niUSRP_Error_UserPayload)

#define niUSRP_Error_Description_MemoryAllocation \
   "Unable to allocate the requested memory. "
#define niUSRP_Error_Description_TimeoutExceeded \
   "Timeout exceeded before packet received or sent.  "\
   "Not all samples may have been received or sent.  "\
   "Consider increasing timeout."
#define niUSRP_Error_Description_LateCommand \
   "A stream command was issued in the past."
#define niUSRP_Error_Description_BrokenChain \
   "Expected another stream command."
#define niUSRP_Error_Description_Overflow \
   "Overflow: an internal receive buffer has filled before the "\
   "data could be returned.  Consider reducing the IQ rate, "\
   "increasing the Fetch rate, or increasing "\
   "the number of samples per Fetch."
#define niUSRP_Error_Description_BadPacket \
   "The packet could not be parsed."
#define niUSRP_Error_Description_Underflow \
   "Underflow: the Tx buffer was emptied before new data "\
   "was provided.  Consider reducing the IQ rate, "\
   "increasing the Write rate, or increasing "\
   "the number of samples per Write."
#define niUSRP_Error_Description_SeqError \
   "Packet loss between host and device. This may be a symptom of overflows due to an inability to "\
   "maintain streaming at the specified IQ Rate, or rearrangement of packets by an ethernet adapter "\
   "or network switch."
#define niUSRP_Error_Description_TimeError \
   "Packet had timestamp that was late (or too early)."
#define niUSRP_Error_Description_PacketUnderflow \
   "Underflow: occurred inside a packet."
#define niUSRP_Error_Description_BurstSeqError \
   "Packet loss within a burst."
#define niUSRP_Error_Description_InvalidSessionHandle \
   "The specified session handle is invalid or does not"\
   " correspond to an active session."
#define niUSRP_Error_Description_InvalidAttributeID \
   "The specified attribute ID is not recognized."\
   " The attribute may not be valid for this session type (Tx vs. Rx)"\
   " or may not be supported by this device model."
#define niUSRP_Error_Description_InvalidAttributeScope \
   "The specified attribute ID is not valid for the"\
   " specified scope (or channel)."
#define niUSRP_Error_Description_UnsupportedDevice \
   "The specified device (or configuration of multiple"\
   " devices) is not supported."
#define niUSRP_Error_Description_InvalidTimeSpecifier \
   "The given time specifier is not recognzied."
#define niUSRP_Error_Description_InvalidTriggerType \
   "The specified trigger type is not valid."
#define niUSRP_Error_Description_InvalidSourceTerminal \
   "The specified source or output terminal is not valid for this attribute."
#define niUSRP_Error_Description_InvalidPolarity \
   "The specified polarity is invalid or not supported for this source attribute."
#define niUSRP_Error_Description_DriverInRunningState \
   "This attribute cannot be modified while the driver is in the Running state."
#define niUSRP_Error_Description_DriverNotInRunningState \
   "This operation requires the driver to be in the Running state."
#define niUSRP_Error_Description_MultipleAttributeGetConflict \
   "The attribute value differs across channels."
#define niUSRP_Error_Description_InvalidChannelName \
   "The specified channel name is not recognized."
#define niUSRP_Error_Description_InvalidDeviceName \
   "The specified device name is not recognized."
#define niUSRP_Error_Description_InvalidFetchSize \
   "The specified number of samples to fetch is invalid for this configuration."
#define niUSRP_Error_Description_InternalDriverError \
   "A runtime or configuration error occurred."
#define niUSRP_Error_Description_LocalOscillatorFailedToLock \
   "The local oscillator did not lock within the allotted time."
#define niUSRP_Error_Description_ReferenceClockFailedToLock \
   "The reference clock PLL did not lock within the allotted time."
#define niUSRP_Error_Description_MimoLinkFailedToLock \
   "The MIMO cable link did not lock within the allotted time."
#define niUSRP_Error_Description_InvalidTimeout \
   "Timeout must be -1 or >= 0."
#define niUSRP_Error_Description_InvalidSampleWidth \
   "Sample width can be 8 or 16."
#define niUSRP_Error_Description_InvalidHostDataType \
   "Host data type can be niUSRP_Val_HostDataType_ComplexDouble"\
   " or niUSRP_Val_HostDataType_ComplexInt16."
#define niUSRP_Error_Description_HostDataTypeConflict \
   "The type of Fetch or Write method invoked does not match the "\
   " niUSRP_Attr_HostDataType value."
#define niUSRP_Error_Description_InvalidSizeOfDataBuffer \
   "The number of channels in the session does not correspond to the "\
   " dimensionality of the data buffer.  Use the 2D Fetch or Write VIs "\
   " when using multi-device sessions."
#define niUSRP_Error_Description_Alignment \
   "The incoming data packets could not be time-aligned.  This could be"\
   " due to devices not being synchronized or failure to stream fast enough"\
   " for specified IQ rate."
#define niUSRP_Error_Description_PpsWithGPSDC \
   "Cannot use the PPS In port while GPSDC attached"
#define niUSRP_Error_Description_GPSFailedToLock \
   "The GPS did not lock within the allotted time. The first time may take up to 10 minutes."
#define niUSRP_Error_Description_UnsupportedLockSensor \
   "The desired sensor lock is unsupported"
#define niUSRP_Error_Description_FunctionNotSupported \
   "This method is not supported by this particular session type (Tx vs. Rx)"\
   " or is not supported by this device model."
#define niUSRP_Error_Description_InvalidOptionString \
   "The option string is improperly formatted or attempts to set unsupported options or values."
#define niUSRP_Error_Description_InvalidSignal \
   "The specified signal is unrecognized or not supported on this device."
#define niUSRP_Error_Description_IQRateDiffersAmongChannels \
   "The IQ Rate must be set to the same value across all channels."
#define niUSRP_Error_Description_OutOfSequence \
   "The received data is out of sequence. This may be a symptom of overflows due to an inability to "\
   "maintain streaming at the specified IQ Rate, or rearrangement of packets by an ethernet adapter "\
   "or network switch."
#define niUSRP_Error_Description_UserPayload \
   "User payload streaming error."
#define niUSRP_Error_Description_Invalid290xDeviceName \
   "The specified device name is not recognized. If using an NI 290x, ensure that only a single device is specified."

#endif /* ___niUSRPErrors_h___ */
