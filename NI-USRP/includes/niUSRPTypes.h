/*
 * Types for the NI-USRP C API
 *
 * Copyright (c) 2011-
 * National Instruments Corporation.
 * All rights reserved.
 */

#ifndef ___niUSRPTypes_h___
#define ___niUSRPTypes_h___

#include <time.h>

/*
 * Determine platform details.
 */
#if defined(_M_IX86) \
 || defined(_M_X64) \
 || defined(i386) \
 || defined(__i386__) \
 || defined(__amd64__) \
 || defined(__amd64) \
 || defined(__x86_64__) \
 || defined(__x86_64) \
 || defined(__i386) \
 || defined(_X86_) \
 || defined(__THW_INTEL__) \
 || defined(__I86__) \
 || defined(__INTEL__) \
 || defined(__X86__) \
 || defined(__386__) \
 || defined(__I86__) \
 || defined(M_I386) \
 || defined(M_I86) \
 || defined(_M_I386) \
 || defined(_M_I86)
   #if defined(_WIN32) \
    || defined(_WIN64) \
    || defined(__WIN32__) \
    || defined(__TOS_WIN__) \
    || defined(__WINDOWS__) \
    || defined(_WINDOWS) \
    || defined(__WINDOWS_386__) \
    || defined(__CYGWIN__)
      /* Either Windows or Phar Lap ETS. */
      #define niUSRP_Windows 1
   #elif defined(__linux) \
      || defined(__linux__) \
      || defined(__gnu_linux__) \
      || defined(linux)
      #define niUSRP_Linux 1
   #else
      #error Unsupported OS.
   #endif
#elif defined(__powerpc) \
   || defined(__powerpc__) \
   || defined(__POWERPC__) \
   || defined(__ppc__) \
   || defined(__PPC) \
   || defined(_M_PPC) \
   || defined(_ARCH_PPC) \
   || defined(__PPC__) \
   || defined(__ppc)
   #if defined(__vxworks)
      #define niUSRP_VxWorks 1
   #else
      #error Unsupported OS.
   #endif
#else
   #error Unsupported architecture.
#endif

/*
 * Determine compiler.
 */
#if defined(_MSC_VER)
   #define niUSRP_Msvc 1
#elif defined(__GNUC__)
   #define niUSRP_Gcc 1
#elif defined(_CVI_) && !defined(_TPC_)
   #define niUSRP_Cvi 1
   /* Enables CVI Library Protection Errors. */
   #pragma EnableLibraryRuntimeChecking
#else
   /* Unknown compiler. */
#endif

/*
 * Determine compliance with different C/C++ language standards.
 */
#if defined(__cplusplus)
   #define niUSRP_Cpp 1
   #if __cplusplus >= 199707L
      #define niUSRP_Cpp98 1
   #endif
#endif
#if defined(__STDC__)
   #define niUSRP_C89 1
   #if defined(__STDC_VERSION__)
      #define niUSRP_C90 1
      #if __STDC_VERSION__ >= 199409L
         #define niUSRP_C94 1
         #if __STDC_VERSION__ >= 199901L
            #define niUSRP_C99 1
         #endif
      #endif
   #endif
#endif

/*
 * Determine ability to inline functions.
 */
#if niUSRP_Cpp || niUSRP_C99
   /* The inline keyword exists in C++ and C99. */
   #define niUSRP_Inline inline
#elif niUSRP_Msvc
   /* Visual C++ (at least since 6.0) also supports an alternate keyword. */
   #define niUSRP_Inline __inline
#elif niUSRP_Gcc
   /* GCC (at least since 2.95.2) also supports an alternate keyword. */
   #define niUSRP_Inline __inline__
#elif !defined(niUSRP_Inline)
   /*
    * Disable inlining if inline support is unknown. To manually enable
    * inlining, #define the following macro before #including niUSRP.h:
    *
    *    #define niUSRP_Inline inline
    */
   #define niUSRP_Inline
#endif

/*
 * Define exact-width integer types, if they have not already been defined.
 */
#if niUSRP_ExactWidthIntegerTypesDefined \
 || defined(_STDINT) \
 || defined(_STDINT_H) \
 || defined(_STDINT_H_) \
 || defined(_INTTYPES_H) \
 || defined(_INTTYPES_H_) \
 || defined(_SYS_STDINT_H) \
 || defined(_SYS_STDINT_H_) \
 || defined(_SYS_INTTYPES_H) \
 || defined(_SYS_INTTYPES_H_) \
 || defined(_STDINT_H_INCLUDED) \
 || defined(BOOST_CSTDINT_HPP) \
 || defined(_MSC_STDINT_H_) \
 || defined(_PSTDINT_H_INCLUDED)
   /* Assume that exact-width integer types have already been defined. */
#elif niUSRP_C99 \
   || niUSRP_Gcc /* GCC (at least since 3.0) has a stdint.h. */ \
   || defined(HAVE_STDINT_H)
   /* Assume that stdint.h can be included. */
   #include <stdint.h>
#elif niUSRP_Msvc \
   || niUSRP_Cvi
   /* Manually define exact-width integer types. */
   typedef   signed    char  int8_t;
   typedef unsigned    char uint8_t;
   typedef            short  int16_t;
   typedef unsigned   short uint16_t;
   typedef              int  int32_t;
   typedef unsigned     int uint32_t;
   typedef          __int64  int64_t;
   typedef unsigned __int64 uint64_t;
#else
   /*
    * Exact-width integer types must be defined by the user, and the following
    * macro must be #defined, before #including niUSRP.h:
    *
    *    #define niUSRP_ExactWidthIntegerTypesDefined 1
    */
   #error Exact-width integer types must be defined by the user. See comment.
#endif

/* Included for definition of size_t. */
#include <stddef.h>

#if niUSRP_Cpp
extern "C" {
#endif

/**
 * String types
 */
typedef char* niUSRP_String;
typedef const char* niUSRP_ConstString;

/**
 * A boolean value; either niUSRP_False or niUSRP_True.
 */
typedef uint16_t niUSRP_Bool;

/**
 * Represents a false condition.
 */
static const niUSRP_Bool niUSRP_False = 0;

/**
 * Represents a true condition.
 */
static const niUSRP_Bool niUSRP_True = 1;

// maximum string length constants for niUSRP_Device struct
static const size_t niUSRP_Val_MaximumAddressLength = 16;
static const size_t niUSRP_Val_MaximumSerialNumberLength = 16;
static const size_t niUSRP_Val_MaximumModelLength = 36;

/**
 * Represents the resulting status of a function call through its return value.
 * 0 is success, negative values are errors, and positive values are warnings.
 */
typedef int32_t niUSRP_Status;
#define mAssignStatus(target,source) \
   target=(target<0)?target:((source<0)?source:target);


/**
 * Session type
 */
typedef uint32_t niUSRP_Session;
static const niUSRP_Session niUSRP_InvalidSession = 0;

/**
 * Attribute ID type
 */
typedef uint32_t niUSRP_AttributeID;

/**
 * Represents a complex double-precision floating point value.
 */
#pragma pack(push,8)
typedef struct niUSRP_ComplexDouble_struct
{
   double real;
   double imaginary;
} niUSRP_ComplexDouble;

typedef struct niUSRP_ComplexInt16_struct
{
   int16_t real;
   int16_t imaginary;
} niUSRP_ComplexInt16;

typedef struct niUSRP_Timestamp_struct
{
   time_t wholeSeconds;
   double fractionalSeconds;
} niUSRP_Timestamp;

typedef struct niUSRP_Device_struct
{
   char address[niUSRP_Val_MaximumAddressLength];
   char serialNumber[niUSRP_Val_MaximumSerialNumberLength];
   char model[niUSRP_Val_MaximumModelLength];
} niUSRP_Device;
#pragma pack(pop)

/**
 * API Export Decoration
 */
#define niUSRP_API_Function __declspec(dllexport) niUSRP_Status __stdcall 
#define niUSRP_API_Function_Internal niUSRP_Status __stdcall 

#if niUSRP_Cpp
}
#endif

#endif /* ___niUSRPTypes_h___ */
