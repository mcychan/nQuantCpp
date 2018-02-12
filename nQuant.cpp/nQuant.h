
// nQuant.cpp.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CQuantApp:
// See nQuant.cpp.cpp for the implementation of this class
//

class CQuantApp : public CWinApp
{

private:
	GdiplusStartupInput  m_gdiplusStartupInput;
	ULONG_PTR m_gdiplusToken;

public:
	CQuantApp();

// Overrides
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CQuantApp theApp;