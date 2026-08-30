// stdafx.h
#pragma once
#define _USE_MATH_DEFINES

#ifndef GDIPVER
	#define GDIPVER 0x0110 // Enable GDI+ 1.1 features
#endif

#define COM_NO_WINDOWS_H

#if defined(_WIN32) || defined(_WIN64)
	#include <unknwn.h>         // Base COM interfaces (IUnknown)
#endif

#if defined(_MSC_VER)
	// --- Microsoft Visual Studio ---
	#pragma comment(lib, "gdiplus")

	#include <gdiplus.h>

	#ifdef _WIN64
		#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='amd64' publicKeyToken='6595b64144ccf1df' language='*'\"")
	#else
		#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='X86' publicKeyToken='6595b64144ccf1df' language='*'\"")
	#endif

#elif defined(__MINGW32__) || defined(__MINGW64__)
	// --- MinGW-w64 (GCC for Windows) ---
	#include <windows.h>
	#include <gdiplus.h>
	#include <objbase.h>

#else
	// --- Wine / WineG++ Compiler
	#include <windows.h>
	#include <guiddef.h>
	#include <gdiplus/gdiplus.h>

#endif

// Shared namespace usage across all environments
using namespace Gdiplus;

#ifndef BYTE_MAX
	#include <climits>
	#include <cmath>
	#define BYTE_MAX UCHAR_MAX
#endif

inline double sqr(double value)
{
	return value * value;
}

#ifdef _WIN64
#define _sqrt sqrt
#elif defined(_WIN32)
inline double __declspec (naked) __fastcall _sqrt(double n)
{
	_asm fld qword ptr[esp + 4]
		_asm fsqrt
	_asm ret 8
}
#endif