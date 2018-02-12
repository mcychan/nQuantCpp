// stdafx.h

#pragma once

#include "targetver.h"

#include <tchar.h>
#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS
#define _AFX_NO_MFC_CONTROLS_IN_DIALOGS

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif

#include <afx.h>
#include <afxwin.h>
#include <afxext.h>
#ifndef _AFX_NO_OLE_SUPPORT
#include <afxdtctl.h>
#endif
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>
#endif // _AFX_NO_AFXCMN_SUPPORT

#include <iostream>


#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment(linker, "/manifestdependency:\"type='win32' name='Microsoft.Windows.GdiPlus' version='1.1.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

