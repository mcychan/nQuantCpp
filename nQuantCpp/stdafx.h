// stdafx.h

#pragma once

#include <atlstr.h>
#include <tchar.h>
#include <intsafe.h>

#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS      // some CString constructors will be explicit

#define GDIPVER 0x0110 //Use more advanced GDI+ features
#pragma comment(lib, "gdiplus.lib")
#include <unknwn.h>				// Needed for non-MFC/ATL use
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

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