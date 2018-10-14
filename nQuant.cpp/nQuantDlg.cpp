
// nQuantDlg.cpp : implementation file
//

#include "stdafx.h"
#include "nQuant.h"
#include "nQuantDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About
class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CQuantDlg dialog
CQuantDlg::CQuantDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_NQUANTCPP_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CQuantDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CQuantDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDRETRY, &CQuantDlg::OnBnClickedRetry)
	ON_BN_CLICKED(ID_FILE_OPEN, &CQuantDlg::OnBnClickedFileOpen)
END_MESSAGE_MAP()


// CQuantDlg message handlers

BOOL CQuantDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CQuantDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CQuantDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CQuantDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

LPCTSTR CQuantDlg::GetInitialDir()
{
	SHGetSpecialFolderPath(NULL, szInitDir, CSIDL_MYPICTURES, FALSE);
	return szInitDir;
}

CString CQuantDlg::GetFileType()
{
	CString strFilter;
	CString strAllSupported = _T("All Supported Formats|");

	int nIndex = 0;
	UINT  nDecCount = 0;  // number of image decoders
	UINT  nDecSize = 0;  // size of the image decoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;

	GetImageDecodersSize(&nDecCount, &nDecSize);
	pImageCodecInfo = (ImageCodecInfo*)malloc(nDecSize);
	GetImageDecoders(nDecCount, nDecSize, pImageCodecInfo);

	for (UINT nIndex = 0; nIndex < nDecCount; nIndex++)
	{
		CString strFilterItem;
		CString fileExtension = CString(pImageCodecInfo[nIndex].FilenameExtension);
		strFilterItem.Format(_T("%s files(%s)|%s|"),
			pImageCodecInfo[nIndex].FormatDescription,
			fileExtension.MakeLower(),
			fileExtension);
		strFilter += strFilterItem;
		strAllSupported += fileExtension;
		if (nIndex < nDecCount - 1)
		{
			strAllSupported += _T(";");
		}
	}
	free(pImageCodecInfo);

	// Prepend "All GDI+ supported" to filter
	m_FileTypes.Format(_T("%s|%s"), strAllSupported, strFilter);

	// Get filter into commdlg format (lots of \0)
	if (!m_FileTypes.IsEmpty())
	{
		LPTSTR pch = m_FileTypes.GetBuffer(0); // modify the buffer in place
											   // MFC delimits with '|' not '\0'
		while ((pch = _tcschr(pch, '|')) != NULL)
			*pch++ = '\0';
	}

	return m_FileTypes;
}

//////////////////
// Load from path name.
//
Bitmap* CQuantDlg::Load(LPCTSTR pszPathName)
{
	IStream* fileStream = NULL;
	HRESULT hr = SHCreateStreamOnFile(
		pszPathName,
		STGM_READ | STGM_SHARE_DENY_NONE,
		&fileStream);

	if (!SUCCEEDED(hr)) {
		return NULL;
	}

	m_PathName = pszPathName;
	Bitmap* result = Bitmap::FromStream(fileStream, TRUE);
	fileStream->Release();
	return result;
}

void CQuantDlg::OnBnClickedFileOpen()
{
	CFileDialog dlgFile(TRUE, NULL, NULL, OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST, NULL, AfxGetMainWnd());
	dlgFile.m_ofn.lpstrInitialDir = GetInitialDir();
	dlgFile.m_ofn.lpstrFilter = GetFileType();

	if (IDOK == dlgFile.DoModal()) {		
		m_pImage.reset(Load(dlgFile.GetPathName()));
		HBITMAP bitmap;
		m_pImage->GetHBITMAP(NULL, &bitmap);
		SetBackgroundImage(bitmap, CDialogEx::BACKGR_TOPLEFT);
	}
}

void CQuantDlg::OnBnClickedRetry()
{
	theApp.DoWaitCursor(1); // 1->>display the hourglass cursor
	UINT w = m_pImage->GetWidth();
	UINT h = m_pImage->GetHeight();

	m_pTargetImage = make_unique<Bitmap>(w, h, PixelFormat8bppIndexed);

	UINT nMaxColors = 256;
	bool bSucceeded = pnnLABQuantizer.QuantizeImage(m_pImage.get(), m_pTargetImage.get(), nMaxColors);
	//bool bSucceeded = easQuantizer.QuantizeImage(m_pImage.get(), m_pTargetImage.get(), nMaxColors);

	HBITMAP bitmap;
	m_pTargetImage->GetHBITMAP(NULL, &bitmap);
	
	TCHAR path[MAX_PATH];
	GetTempPath(MAX_PATH, path);
	CString pathName;
	pathName.Format(_T("%s%s_new.png"), path, PathFindFileName(m_PathName));

	// image/png  : {557cf406-1a04-11d3-9a73-0000f81ef32e}
	const CLSID pngEncoderClsId = { 0x557cf406, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
	Status status = m_pTargetImage->Save(pathName, &pngEncoderClsId);
	if(bSucceeded && status != Status::Ok)
		SetBackgroundImage(bitmap, CDialogEx::BACKGR_TOPLEFT);
	theApp.DoWaitCursor(-1); // -1->>remove the hourglass cursor
	if (bSucceeded) {
		CString msg = _T("Converted to 256 colors successfully");
		if (status == Status::Ok) {
			msg += _T(" [") + pathName + _T("]");
			/*HINSTANCE hInst = */ShellExecute(GetSafeHwnd(), _T("open"), pathName, NULL, NULL, SW_SHOWNORMAL);
		}
		AfxMessageBox(msg, MB_ICONINFORMATION);
	}
	else
		AfxMessageBox(_T("Failed to convert to 256 colors"));
}
