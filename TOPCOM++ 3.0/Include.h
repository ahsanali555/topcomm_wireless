/*! \file 
\brief platform dependent type definitions
*/
#pragma once
#ifdef WIN32
	#include <conio.h>
	typedef unsigned __int64 uint64;
	#define forceinline __forceinline
	#define accessg _access
#else
	#include <sys/time.h>
	#include <stdio.h>
	typedef unsigned long long  uint64;
	typedef long long __int64;
	int _kbhit(void);
	int _getch(void);
	#define forceinline inline
	#define accessg access
#endif

