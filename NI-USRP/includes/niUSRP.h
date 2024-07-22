/*****************************************************************************
 *===========================================================================*
 *    Copyright 2011, National Instruments, Corporation. All Rights Reserved.*
 *****************************************************************************/

/****************************************************************************
 *                       niUSRP
 *
 * Title:    niUSRP.h
 * Purpose:  niUSRP
 *           instrument driver declarations.
 *
 ****************************************************************************/

#ifndef ___niUSRP_h___
#define ___niUSRP_h___

#ifndef ___niUSRPTypes_h___
   #include "niUSRPTypes.h"
#endif
#ifndef ___niUSRPErrors_h___
   #include "niUSRPErrors.h"
#endif

#if niUSRP_Cpp
extern "C" {
#endif


/****************************************************************************
 *  Attributes
 ****************************************************************************/

#define niUSRP_Attr_Base                                       1150000

#define niUSRP_Attr_NumberOfSamplesIsFinite                    (niUSRP_Attr_Base+0x01) /* niUSRP_Bool */
#define niUSRP_Attr_NumberOfSamples                            (niUSRP_Attr_Base+0x02) /* int64_t */
#define niUSRP_Attr_IQRate                                     (niUSRP_Attr_Base+0x03) /* double:    Samples/s */
#define niUSRP_Attr_CarrierFrequency                           (niUSRP_Attr_Base+0x04) /* double:    Hz */
#define niUSRP_Attr_Bandwidth                                  (niUSRP_Attr_Base+0x05) /* double:    Hz */
#define niUSRP_Attr_Gain                                       (niUSRP_Attr_Base+0x06) /* double */
#define niUSRP_Attr_ActiveAntenna                              (niUSRP_Attr_Base+0x07) /* niUSRP_String: enum niUSRP_Terminal */
#define niUSRP_Attr_ReferenceFrequencySource                   (niUSRP_Attr_Base+0x08) /* niUSRP_String: enum niUSRP_Terminal */
#define niUSRP_Attr_TimebaseClockSource                        (niUSRP_Attr_Base+0x09) /* niUSRP_String: enum niUSRP_Terminal */
#define niUSRP_Attr_TimebaseClockPolarity                      (niUSRP_Attr_Base+0x0A) /* int32_t:   enum niUSRP_SignalPolarity */
#define niUSRP_Attr_DriverState                                (niUSRP_Attr_Base+0x0B) /* int32_t:   enum niUSRP_DriverState */
#define niUSRP_Attr_SampleWidth                                (niUSRP_Attr_Base+0x0C) /* int32_t */
#define niUSRP_Attr_ExpectedPeak                               (niUSRP_Attr_Base+0x0D) /* double:    fraction of full-scale */
#define niUSRP_Attr_LOFrequency                                (niUSRP_Attr_Base+0x0E) /* double:    Hz */
#define niUSRP_Attr_HostDataType                               (niUSRP_Attr_Base+0x0F) /* int32_t:   enum niUSRP_HostDataType */
#define niUSRP_Attr_EnabledChannels                            (niUSRP_Attr_Base+0x10) /* niUSRP_String */

#define niUSRP_Attr_StartTriggerType                           (niUSRP_Attr_Base+0x50) /* int32_t:   enum niUSRP_TriggerType */
#define niUSRP_Attr_StartTriggerTimeWhole                      (niUSRP_Attr_Base+0x51) /* int64_t */
#define niUSRP_Attr_StartTriggerTimeFractional                 (niUSRP_Attr_Base+0x52) /* double */
#define niUSRP_Attr_ExportedReferenceFrequencyOutputTerminal   (niUSRP_Attr_Base+0x53) /* niUSRP_String:   enum niUSRP_Terminal */
#define niUSRP_Attr_ExportedTimebaseClockOutputTerminal        (niUSRP_Attr_Base+0x54) /* niUSRP_String:   enum niUSRP_Terminal */

#define niUSRP_Attr_CurrentDriverVersion                       (niUSRP_Attr_Base+0x100)/* niUSRP_String:   read-only */
#define niUSRP_Attr_Model                                      (niUSRP_Attr_Base+0x101)/* niUSRP_String:   read-only */
#define niUSRP_Attr_CurrentFirmwareVersion                     (niUSRP_Attr_Base+0x102)/* niUSRP_String:   read-only */
#define niUSRP_Attr_OldestCompatibleFirmwareVersion            (niUSRP_Attr_Base+0x103)/* niUSRP_String:   read-only */
#define niUSRP_Attr_CurrentFPGAVersion                         (niUSRP_Attr_Base+0x104)/* niUSRP_String:   read-only */
#define niUSRP_Attr_OldestCompatibleFPGAVersion                (niUSRP_Attr_Base+0x105)/* niUSRP_String:   read-only */
#define niUSRP_Attr_HasGPSDC                                   (niUSRP_Attr_Base+0x106)/* niUSRP_Bool:   read-only */

#define niUSRP_Attr_GPSTime                                    (niUSRP_Attr_Base+0x130)/* int64_t:   time in seconds since Jan 1, 1970 */
#define niUSRP_Attr_GPSLockStatus                              (niUSRP_Attr_Base+0x131)/* niUSRP_Bool */
#define niUSRP_Attr_GPSSentenceGGA                             (niUSRP_Attr_Base+0x132)/* niUSRP_String */
#define niUSRP_Attr_GPSSentenceRMC                             (niUSRP_Attr_Base+0x133)/* niUSRP_String */

#define niUSRP_Attr_NumberOfChannelsInSession                  (niUSRP_Attr_Base+0x150)/* int32_t:   read-only */
#define niUSRP_Attr_NumberOfMotherboardsInSession              (niUSRP_Attr_Base+0x151)/* int32_t:   read-only */

#define niUSRP_Attr_WarningPolicy                              (niUSRP_Attr_Base+0x200)/* int32_t:   enum niUSRP_WarningPolicy */

/****************************************************************************
 *  Enumerations
 ****************************************************************************/

/**
 * enum niUSRP_TriggerType
 * Used by attribute niUSRP_Attr_StartTriggerType
 */
#define niUSRP_Val_TriggerType_None                   0 /* No trigger (Immediate) */
#define niUSRP_Val_TriggerType_Time                   1 /* Synchronous to a timestamp */

/**
 * enum niUSRP_SignalPolarity
 * Used by attribute niUSRP_Attr_TimebaseClockPolarity
 */
#define niUSRP_Val_SignalPolarity_Rising              0 /* Rising Edge */
#define niUSRP_Val_SignalPolarity_Falling             1 /* Falling Edge */

/**
 * enum niUSRP_Time
 * Used by the functions niUSRP_SetTime and niUSRP_GetTime
 */
#define niUSRP_Val_Time_Now                           0 /* Now */
#define niUSRP_Val_Time_NextTimebaseEdge              1 /* At the next timebase clock edge */
#define niUSRP_Val_Time_PreviousTimebaseEdge          2 /* The previous timebase clock edge */

/**
 * enum niUSRP_Terminal
 * Used by the attributes niUSRP_ReferenceFrequencySource and niUSRP_Attr_TimebaseClockSource
 */
#define niUSRP_Val_Terminal_Internal         "Internal"  /* The internal reference */
#define niUSRP_Val_Terminal_RefIn            "RefIn"     /* The REF IN front-panel terminal */
#define niUSRP_Val_Terminal_PpsIn            "PpsIn"     /* The PPS IN front-panel terminal */
#define niUSRP_Val_Terminal_Mimo             "Mimo"      /* The MIMO EXPANSION front-panel terminal */
#define niUSRP_Val_Terminal_Rx1              "RX1"       /* The RX1 front-panel terminal */
#define niUSRP_Val_Terminal_Rx2              "RX2"       /* The RX2 front-panel terminal */
#define niUSRP_Val_Terminal_Tx1              "TX1"       /* The RX1 front-panel terminal */
#define niUSRP_Val_Terminal_Tx2              "TX2"       /* The RX2 front-panel terminal */
#define niUSRP_Val_Terminal_GPS              "GPS"       /* A GPS is attached to the external reference */
#define niUSRP_Val_Terminal_RefOut           "RefOut"    /* The REF OUT front-panel terminal */
#define niUSRP_Val_Terminal_PpsTrigOut       "PpsTrigOut"/* The PPS TRIG OUT front-panel terminal */

/**
 * enum niUSRP_Signal
 * Used by the function niUSRP_ExportSignal
 */
#define niUSRP_Val_Signal_ReferenceFrequency    0     /* The reference frequency signal */
#define niUSRP_Val_Signal_TimebaseClock         1     /* The timebase (often PPS) signal */

/**
 * enum niUSRP_DriverState
 * Used by the attribute niUSRP_Attr_DriverState
 */
#define niUSRP_Val_DriverState_Initialized            1 /* Initialized */
#define niUSRP_Val_DriverState_Committed              4 /* Committed */
#define niUSRP_Val_DriverState_Running                5 /* Running */
#define niUSRP_Val_DriverState_Done                   6 /* Done */

/**
 * enum niUSRP_HostDataType
 * Used by the attribute niUSRP_Attr_HostDataType
 */
#define niUSRP_Val_HostDataType_ComplexDouble            0 /* niUSRP_ComplexDouble */
#define niUSRP_Val_HostDataType_ComplexInt16             1 /* niUSRP_ComplexInt16 */

/**
 * enum niUSRP_WarningPolicy
 * Used by the attribute niUSRP_Attr_WarningPolicy
 */
#define niUSRP_Val_WarningPolicy_ReturnWarnings       0 /* Warnings can be returned */
#define niUSRP_Val_WarningPolicy_ReturnErrors         1 /* Errors substituted for warnings */
#define niUSRP_Val_WarningPolicy_ReturnSuccess        2 /* Success substituted for warnings */


/****************************************************************************
 *  Functions
 ****************************************************************************/

/**
 * Creates a new Rx session to the device(s).
 * @param deviceNames Name(s) or IP address(es) of the device(s).
 * @param reset Specifies whether or not to reset the device(s) to a known initialization state.
 * @param sessionHandle Newly created session handle used in subsequent driver calls.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_OpenRxSession(
   niUSRP_String deviceNames,
   niUSRP_Bool reset,
   niUSRP_Session* sessionHandle);

/**
 * Creates a new Tx session to the device(s).
 * @param deviceNames Name(s) or IP address(es) of the device(s).
 * @param reset Specifies whether or not to reset the device(s) to a known initialization state.
 * @param sessionHandle Newly created session handle used in subsequent driver calls.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_OpenTxSession(
   niUSRP_String deviceNames,
   niUSRP_Bool reset,
   niUSRP_Session* sessionHandle);

/**
 * Creates a new session with options to the device(s).
 * @param deviceNames Name(s) or IP address(es) of the device(s).
 * @param options Specifies session settings.
 * @param reset Specifies whether or not to reset the device(s) to a known initialization state.
 * @param sessionHandle Newly created session handle used in subsequent driver calls.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_OpenSessionWithOptions(
   niUSRP_String deviceNames,
   niUSRP_String options,
   niUSRP_Bool reset,
   niUSRP_Session* sessionHandle);


/**
 * Closes the session to the device(s).
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_CloseSession(
   niUSRP_Session sessionHandle);

/**
 * Configures several basic properties of the Tx or Rx signal.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) to which to apply the new settings.
 * @param iqRate Specifies the rate (Samples/s) of the baseband IQ data.
 * @param carrierFrequency Specifies the center frequency (Hz) of the RF signal.
 * @param gain Specifies the aggregate gain applied to the RF signal.
 * @param activeAntenna Specifies the antenna to be used.
 * @param coercedIqRate The actual IQ rate used for this session (coerced to the capabilities of the device).
 * @param coercedCarrierFrequency The actual carrier frequency used for this session (coerced to the capabilities of the device).
 * @param coercedGain The actual gain used for this session (coerced to the capabilities of the device).
 * @return result of the call
 */
niUSRP_API_Function niUSRP_ConfigureSignal(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   double iqRate,
   double carrierFrequency,
   double gain,
   niUSRP_ConstString activeAntenna,
   double* coercedIqRate,
   double* coercedCarrierFrequency,
   double* coercedGain);

/**
 * Configures the start trigger to be generated from the onboard timer and specifies the time the trigger should occur.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param startTriggerTime Specifies the time that the trigger should occur.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_ConfigureTimeStartTrigger(
   niUSRP_Session sessionHandle,
   niUSRP_Timestamp startTriggerTime);

/**
 * Disables the start trigger. When the start trigger is disabled, acquisitions and generations will begin immediately when data is available.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_DisableStartTrigger(
   niUSRP_Session sessionHandle);

/**
 * Sets the time value of the onboard timer.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) to which to apply the new settings.
 * @param when Specifies when to set the specified timestamp on the onboard timer.
 * @param timestamp Specifies the time value to set on the onboard timer.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetTime(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int32_t when,
   niUSRP_Timestamp timestamp);

/**
 * Gets the time value of the onboard timer.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel to which to read the time from
 * @param when Specifies if the timestamp is the time "now" or the time of the last pps.
 * @param timestamp Indicates the fetched time of the onboard timer.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetTime(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int32_t when,
   niUSRP_Timestamp* timestamp);

/**
 * Tells the device to execute the next hardware configuration setting at the specified timestamp
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to do this on
 * @param timestamp when the command should be run
 */
niUSRP_API_Function niUSRP_SetCommandTime(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_Timestamp timestamp);

/**
 * Tells the device to execute hardware configuration settings immediately
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to do this on
 */
niUSRP_API_Function niUSRP_ClearCommandTime(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList);

/**
 * Configures whether the device operation is finite or continuous and the number of samples.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param numberOfSamplesIsFinite Specifies whether the device will generate/acquire a fixed number of samples or generate/acquire continuously.
 * @param numberOfSamples Specifies the number of samples that will be generated from or acquired by the device.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_ConfigureNumberOfSamples(
						   niUSRP_Session sessionHandle,
						   niUSRP_Bool numberOfSamplesIsFinite,
						   int64_t numberOfSamples);

/**
 * Starts the Rx acquisition. If this method is successful, the driver will be left in the Running state. (This call will also implicitly transition through the Committed state, if necessary). If a Start Trigger has been configured, the device will begin waiting for the specified trigger. 
 Use this function in conjunction with one of the Fetch functions to retrieve the acquired data.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_Initiate(
   niUSRP_Session sessionHandle);

/**
 * Stops an acquisition previously started with niUSRP_Initiate. Unless you want to stop an acquistion before it is complete, calling this function is optional. If this method is successful, the driver will be left in the Done state.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_Abort(
   niUSRP_Session sessionHandle);

/**
 * Validates the attribute values and applies the configuration to hardware.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_Commit(
   niUSRP_Session sessionHandle);

/**
 * Fetches ComplexDouble (complex double precision floating point) data from the specified channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) from which to fetch the data.
 * @param numberOfSamples The number of samples to fetch from the acquisition channel.
 * @param timeout Specifies the time to wait before returning an error if the requested number of sample have not been acquired. A negative value indicates to the driver to wait indefinitely.
 * @param dataBuffer The data buffer that will receive the Rx samples.
 * @param numberOfSamplesReturned The number of data samples actually returned in the data bufffer.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_FetchRxDataComplexDouble(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int64_t numberOfSamples,
   double timeout,
   niUSRP_ComplexDouble dataBuffer[],
   int64_t* numberOfSamplesReturned,
   niUSRP_Timestamp* timestamp);

/**
 * Fetches ComplexInt16 (complex 16-bit integer) data from the specified channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) from which to fetch the data.
 * @param numberOfSamples The number of samples to fetch from the acquisition channel.
 * @param timeout Specifies the time to wait before returning an error if the requested number of sample have not been acquired. A negative value indicates to the driver to wait indefinitely.
 * @param dataBuffer The data buffer that will receive the Rx samples.
 * @param numberOfSamplesReturned The number of data samples actually returned in the data bufffer.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_FetchRxDataComplexInt16(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int64_t numberOfSamples,
   double timeout,
   niUSRP_ComplexInt16 dataBuffer[],
   int64_t* numberOfSamplesReturned,
   niUSRP_Timestamp* timestamp);

/**
 * Writes ComplexDouble (complex double precision floating point) to the specified channel(s).
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) to which to write the data.
 * @param numberOfSamples The number of samples to write to the Tx channel.
 * @param timeout Specifies the time to wait before returning an error if the requested number of sample have not been generated. A negative value indicates to the driver to wait indefinitely.
 * @param dataBuffer The data buffer that holds the Tx samples.
 * @param endOfData Specifies whether this is the last of the data that should be generated, or if more data will be provided for the current generation.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_WriteTxDataComplexDouble(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int64_t numberOfSamples,
   double timeout,
   niUSRP_ComplexDouble dataBuffer[],
   niUSRP_Bool endOfData);

/**
 * Writes ComplexInt16 (complex 16-bit integer) to the specified channel(s).
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) to which to write the data.
 * @param numberOfSamples The number of samples to write to the Tx channel.
 * @param timeout Specifies the time to wait before returning an error if the requested number of sample have not been generated. A negative value indicates to the driver to wait indefinitely.
 * @param dataBuffer The data buffer that holds the Tx samples.
 * @param endOfData Specifies whether this is the last of the data that should be generated, or if more data will be provided for the current generation.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_WriteTxDataComplexInt16(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int64_t numberOfSamples,
   double timeout,
   niUSRP_ComplexInt16 dataBuffer[],
   niUSRP_Bool endOfData);

/**
 * Sets the value of the specified attribute on the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) on which to set the attribute
 * @param attributeID The attribute ID of the attribute to be set
 * @param value The value to assign to the attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetAttributeInt32(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   int32_t value);

/**
 * Sets the value of the specified attribute on the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) on which to set the attribute
 * @param attributeID The attribute ID of the attribute to be set
 * @param value The value to assign to the attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetAttributeInt64(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   int64_t value);

/**
 * Sets the value of the specified attribute on the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) on which to set the attribute
 * @param attributeID The attribute ID of the attribute to be set
 * @param value The value to assign to the attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetAttributeDouble(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   double value);

/**
 * Sets the value of the specified attribute on the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) on which to set the attribute
 * @param attributeID The attribute ID of the attribute to be set
 * @param value The value to assign to the attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetAttributeString(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   niUSRP_ConstString value);

/**
 * Sets the value of the specified attribute on the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) on which to set the attribute
 * @param attributeID The attribute ID of the attribute to be set
 * @param value The value to assign to the attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SetAttributeBool(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   niUSRP_Bool value);

/**
 * Gets the value of the specified attribute from the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param attributeID The ID of the attribute to be queried
 * @param value The value of the specified attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetAttributeInt32(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   int32_t* value);

/**
 * Gets the value of the specified attribute from the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param attributeID The ID of the attribute to be queried
 * @param value The value of the specified attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetAttributeInt64(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   int64_t* value);

/**
 * Gets the value of the specified attribute from the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param attributeID The ID of the attribute to be queried
 * @param value The value of the specified attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetAttributeDouble(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   double* value);

/**
 * Gets the value of the specified attribute from the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param attributeID The ID of the attribute to be queried
 * @param value The value of the specified attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetAttributeString(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   int32_t bufferSize,
   niUSRP_String buffer);

/**
 * Gets the value of the specified attribute from the specified session and channel.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param attributeID The ID of the attribute to be queried
 * @param value The value of the specified attribute
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetAttributeBool(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_AttributeID attributeID,
   niUSRP_Bool* value);

/**
 * Returns the error code and description string for the last error that occurred on the session.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param errorCode The error code for the last error that occurred
 * @param errorDescriptionBufferSize The size of the buffer into which to store the error description
 * @param errorDescriptionBuffer The buffer into which to store the error description
 * @return result of the call
 */
niUSRP_API_Function niUSRP_GetError(
   niUSRP_Session sessionHandle,
   niUSRP_Status* errorCode,
   int32_t errorDescriptionBufferSize,
   niUSRP_String errorDescriptionBuffer);

/**
 * Clears any existing error associated with the session.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_ClearError(
   niUSRP_Session sessionHandle);

/**
 * Resets the device to a known initialization state.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @return result of the call
 */
niUSRP_API_Function niUSRP_Reset(
   niUSRP_Session sessionHandle);

/**
 * Performs a self-test of the device.
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param testResult The result of the self-test.
 * @return result of the call
 */
niUSRP_API_Function niUSRP_SelfTest(
   niUSRP_Session sessionHandle,
   int32_t* testResult);

/**
 * Finds all USRP devices on the system.
 * @param devicesArraySize size of array that will be used as a buffer for the devices
 * @param devices the devices found in the system. If devicesArraySize is >0 but less than the number of devices in the system, this array will store the first devicesArraySize devices found.
 * @param numberOfDevicesInSystem pointer to the number of USRP devices found on the system
 * @return result of the call
 */
niUSRP_API_Function niUSRP_FindDevices(
   int32_t devicesArraySize,
   niUSRP_Device* devices,
   int32_t* numberOfDevicesInSystem);

/**
 * Waits for the specified device/sensor to lock
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channel(s) scope from which to query the attribute
 * @param whatToWaitFor Name of the device/sensor (currently supports: "GPS")
 * @param timeoutSeconds How long to wait in seconds before stopping execution
 * @return result of the call
 */
niUSRP_API_Function niUSRP_WaitForLock(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   niUSRP_ConstString whatToWaitFor,
   double timeout);

/**
 * Checks how many register write packets are waiting in the queue
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param numberOfPkts a pointer passing back the number of packets
 * @return result of the call
 */
niUSRP_API_Function_Internal niUSRP_GetNumberOfOutstandingRegisterPackets(
   niUSRP_Session sessionHandle,
   uint32_t *numberOfPkts);

/**
 * Gets the register writes and formats them into a 2D array
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param packetBufferSize the size of the packet buffer
 * @param packetBuffer a 2D array that will hold the register packets
 * @param numberOfPkts a pointer passing back the number of packets
 * @return result of the call
 */
niUSRP_API_Function_Internal niUSRP_GetRegisterPackets(
   niUSRP_Session sessionHandle,
   uint32_t packetBufferSize,
   uint32_t packetBuffer[],
   uint32_t *numberOfPkts);

/**
 * Creates a route to export the specified signal on the specified output terminal
 * @param sessionHandle Session handle created at niUSRP_OpenRxSession or niUSRP_OpenTxSession
 * @param channelList The channels (in this case, devices) on which to export the signal
 * @param signal The signal to be exported
 * @param signalIdentifier Additional specification of which of multiple signals to export
 * @param outputTerminal The physical terminal on which to output the specified signal
 * @return result of the call
 */
niUSRP_API_Function niUSRP_ExportSignal(
   niUSRP_Session sessionHandle,
   niUSRP_ConstString channelList,
   int32_t signal,
   niUSRP_ConstString signalIdentifier,
   niUSRP_ConstString outputTerminal);

 /****************************************************************************
 *---------------------------- End Include File ----------------------------*
 ****************************************************************************/
#if niUSRP_Cpp
}
#endif

#endif /* ___niUSRP_h___ */
