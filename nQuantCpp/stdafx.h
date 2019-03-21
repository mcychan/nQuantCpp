// stdafx.h

#pragma once

#include "targetver.h"

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
