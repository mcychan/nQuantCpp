// nQuantCpp.cpp
//

#include "stdafx.h"
#include "nQuantCpp.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

GdiplusStartupInput  m_gdiplusStartupInput;
ULONG_PTR m_gdiplusToken;
NeuralNet::NeuQuantizer neuQuantizer;
nQuant::WuQuantizer wuQuantizer;

Bitmap* ConvertTo(Bitmap* pSource, PixelFormat format)
{
	UINT w = pSource->GetWidth();
	UINT h = pSource->GetHeight();

	// If necessary, convert source to 32 bpp RGB bitmap
	if (pSource->GetPixelFormat() == format)
		return pSource;

	auto pSrcBitmap = new Bitmap(w, h, format);

	Graphics g(pSrcBitmap);
	if (g.DrawImage(pSource, 0, 0, w, h) != Ok)
		return NULL;

	return pSrcBitmap;
}

int main(int argc, char** argv)
{
	if (argc < 4 || strcmp(argv[2], "/o")) {
		cout << " nQuantCpp myimage.png /o mynewimage.png";
		return 0;
	}

	VERIFY(GdiplusStartup(&m_gdiplusToken, &m_gdiplusStartupInput, NULL) == Ok);
	{
		unique_ptr<Bitmap> m_pImage;
		unique_ptr<Bitmap> m_pImage256Color;

		TCHAR szDirectory[MAX_PATH];
		::GetCurrentDirectory(sizeof(szDirectory) - 1, szDirectory);
		CString szDir = szDirectory + CString(_T("\\"));
		CString szFile = CString(argv[1]);
		if (szFile.FindOneOf(_T("\\")) < 0)
			szFile = szDir + CString(argv[1]);

		m_pImage.reset(Bitmap::FromFile(szFile));
		Status status = m_pImage->GetLastStatus();
		if (status == Ok) {

			UINT w = m_pImage->GetWidth();
			UINT h = m_pImage->GetHeight();

			// Create 8 bpp indexed bitmap of the same size
			m_pImage256Color = make_unique<Bitmap>(w, h, PixelFormat8bppIndexed);
			Bitmap* pBitmap = ConvertTo(m_pImage.get(), PixelFormat32bppARGB);
			if (pBitmap != m_pImage.get())
				m_pImage.reset(pBitmap);

			UINT nMaxColors = 1 << GetPixelFormatSize(m_pImage256Color->GetPixelFormat());
			bool bSucceeded = neuQuantizer.QuantizeImage(m_pImage.get(), m_pImage256Color.get(), nMaxColors);

			CString pathName = CString(argv[3]);
			if(pathName.FindOneOf(_T("\\")) < 0)
				pathName = szDir + CString(argv[3]);
			// image/png  : {557cf406-1a04-11d3-9a73-0000f81ef32e}
			const CLSID pngEncoderClsId = { 0x557cf406, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
			status = m_pImage256Color->Save(pathName, &pngEncoderClsId);
			if (status != Status::Ok)
				cout << "Failed to save image in '" << argv[3] << "' file";
		}
		else
			cout << "Failed to read image in '" << argv[1] << "' file";
	}
	GdiplusShutdown(m_gdiplusToken);
    return 0;
}
