// nQuantDlg.h : header file
//

#pragma once
#include <memory>
using namespace std;
#include "afxwin.h"
#include "../nQuantCpp/PnnQuantizer.h"
#include "../nQuantCpp/PnnLABQuantizer.h"
#include "../nQuantCpp/EdgeAwareSQuantizer.h"
#include "../nQuantCpp/WuQuantizer.h"

// CQuantDlg dialog
class CQuantDlg : public CDialogEx
{
// Construction
public:
	CQuantDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_NQUANTCPP_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

private:
	LPCTSTR GetInitialDir();
	CString GetFileType();
	Bitmap* Load(LPCTSTR pszPathName);

// Implementation
protected:
	HICON m_hIcon;
	TCHAR szInitDir[MAX_PATH];
	CString m_PathName, m_FileTypes;
	unique_ptr<Bitmap> m_pImage;
	unique_ptr<Bitmap> m_pTargetImage;
	CWinThread* m_lpActionThread;

	PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
	PnnQuant::PnnQuantizer pnnQuantizer;
	nQuant::WuQuantizer wuQuantizer;
	EdgeAwareSQuant::EdgeAwareSQuantizer easQuantizer;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
	DECLARE_MESSAGE_MAP()
public:
	virtual ~CQuantDlg();
	afx_msg void OnBnClickedFileOpen();
	afx_msg void OnBnClickedRetry();
	static UINT DoConvert(LPVOID pParam);
	LRESULT OnFinished(WPARAM wParam, LPARAM lParam);
	
};
