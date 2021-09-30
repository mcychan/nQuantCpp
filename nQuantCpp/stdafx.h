// stdafx.h

#pragma once

#define GDIPVER 0x0110 //Use more advanced GDI+ features
#pragma comment(lib, "gdiplus.lib")
#include <unknwn.h>				// Needed for non-MFC/ATL use
#include <gdiplus.h>
using namespace Gdiplus;
#ifdef _WIN64
	#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='amd64' publicKeyToken='6595b64144ccf1df' language='*'\"")
#else
	#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='X86' publicKeyToken='6595b64144ccf1df' language='*'\"")
#endif

#ifndef BYTE_MAX
	#include <climits>
	#define BYTE_MAX UCHAR_MAX
#endif

#ifndef SHORT_MAX
	#define SHORT_MAX SHRT_MAX
#endif

inline double sqr(double value)
{
	return value * value;
}

#ifdef _WIN64
#define _sqrt sqrt
#else
inline double __declspec (naked) __fastcall _sqrt(double n)
{
	_asm fld qword ptr[esp + 4]
		_asm fsqrt
	_asm ret 8
}
#endif // _WIN64