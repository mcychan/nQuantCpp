#include "stdafx.h"
//////////////////////////////////////////////////////////////////////////
//
// GetBitmapHeaderSize
//
#include "bitmapUtilities.h"

ULONG GetBitmapHeaderSize(LPCVOID pDib)
{
	ULONG nHeaderSize = *(PDWORD)pDib;

	switch (nHeaderSize)
	{
		case sizeof(BITMAPCOREHEADER) :
			case sizeof(BITMAPINFOHEADER) :
			case sizeof(BITMAPV4HEADER) :
			case sizeof(BITMAPV5HEADER) :
		{
			return nHeaderSize;
		}
	}

	return 0;
}


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapLineWidthInBytes
//

ULONG GetBitmapLineWidthInBytes(ULONG nWidthInPixels, ULONG nBitCount)
{
	return (((nWidthInPixels * nBitCount) + 31) & ~31) >> 3;
}


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapDimensions
//

BOOL GetBitmapDimensions(LPCVOID pDib, UINT* pWidth, UINT* pHeight)
{
	ULONG nHeaderSize = GetBitmapHeaderSize(pDib);

	if (nHeaderSize == 0)
	{
		return FALSE;
	}

	if (nHeaderSize == sizeof(BITMAPCOREHEADER))
	{
		PBITMAPCOREHEADER pbmch = (PBITMAPCOREHEADER)pDib;

		if (pWidth != NULL)
		{
			*pWidth = pbmch->bcWidth;
		}

		if (pHeight != NULL)
		{
			*pHeight = pbmch->bcHeight;
		}
	}
	else
	{
		PBITMAPINFOHEADER pbmih = (PBITMAPINFOHEADER)pDib;

		if (pWidth != NULL)
		{
			*pWidth = pbmih->biWidth;
		}

		if (pHeight != NULL)
		{
			*pHeight = abs(pbmih->biHeight);
		}
	}

	return TRUE;
}


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapSize
//

ULONG GetBitmapSize(LPCVOID pDib)
{
	ULONG nHeaderSize = GetBitmapHeaderSize(pDib);

	if (nHeaderSize == 0)
	{
		return 0;
	}

	// Start the calculation with the header size

	ULONG nDibSize = nHeaderSize;

	// is this an old style BITMAPCOREHEADER?

	if (nHeaderSize == sizeof(BITMAPCOREHEADER))
	{
		PBITMAPCOREHEADER pbmch = (PBITMAPCOREHEADER)pDib;

		// Add the color table size

		if (pbmch->bcBitCount <= 8)
		{
			nDibSize += (ULONG)sizeof(RGBTRIPLE) * (1 << pbmch->bcBitCount);
		}

		// Add the bitmap size

		ULONG nWidth = GetBitmapLineWidthInBytes(pbmch->bcWidth, pbmch->bcBitCount);

		nDibSize += nWidth * pbmch->bcHeight;
	}
	else
	{
		// this is at least a BITMAPINFOHEADER

		PBITMAPINFOHEADER pbmih = (PBITMAPINFOHEADER)pDib;

		// Add the color table size

		if (pbmih->biClrUsed != 0)
		{
			nDibSize += sizeof(RGBQUAD) * pbmih->biClrUsed;
		}
		else if (pbmih->biBitCount <= 8)
		{
			nDibSize += (ULONG)sizeof(RGBQUAD) * (1 << pbmih->biBitCount);
		}

		// Add the bitmap size

		if (pbmih->biSizeImage != 0)
		{
			nDibSize += pbmih->biSizeImage;
		}
		else
		{
			// biSizeImage must be specified for compressed bitmaps

			if (pbmih->biCompression != BI_RGB &&
				pbmih->biCompression != BI_BITFIELDS)
			{
				return 0;
			}

			ULONG nWidth = GetBitmapLineWidthInBytes(pbmih->biWidth, pbmih->biBitCount);

			nDibSize += nWidth * abs(pbmih->biHeight);
		}

		// Consider special cases

		if (nHeaderSize == sizeof(BITMAPINFOHEADER))
		{
			// If this is a 16 or 32 bit bitmap and BI_BITFIELDS is used, 
			// bmiColors member contains three DWORD color masks.
			// For V4 or V5 headers, this info is included the header

			if (pbmih->biCompression == BI_BITFIELDS)
			{
				nDibSize += 3 * sizeof(DWORD);
			}
		}
		else if (nHeaderSize >= sizeof(BITMAPV5HEADER))
		{
			// If this is a V5 header and an ICM profile is specified,
			// we need to consider the profile data size

			PBITMAPV5HEADER pbV5h = (PBITMAPV5HEADER)pDib;

			// if there is some padding before the profile data, add it

			if (pbV5h->bV5ProfileData > nDibSize)
			{
				nDibSize = pbV5h->bV5ProfileData;
			}

			// add the profile data size

			nDibSize += pbV5h->bV5ProfileSize;
		}
	}

	return nDibSize;
}


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapOffsetBits
//

ULONG GetBitmapOffsetBits(LPCVOID pDib)
{
	ULONG nHeaderSize = GetBitmapHeaderSize(pDib);

	if (nHeaderSize == 0)
	{
		return 0;
	}

	// Start the calculation with the header size

	ULONG nOffsetBits = nHeaderSize;

	// is this an old style BITMAPCOREHEADER?

	if (nHeaderSize == sizeof(BITMAPCOREHEADER))
	{
		PBITMAPCOREHEADER pbmch = (PBITMAPCOREHEADER)pDib;

		// Add the color table size

		if (pbmch->bcBitCount <= 8)
		{
			nOffsetBits += (ULONG)sizeof(RGBTRIPLE) * (1 << pbmch->bcBitCount);
		}
	}
	else
	{
		// this is at least a BITMAPINFOHEADER

		PBITMAPINFOHEADER pbmih = (PBITMAPINFOHEADER)pDib;

		// Add the color table size

		if (pbmih->biClrUsed != 0)
		{
			nOffsetBits += sizeof(RGBQUAD) * pbmih->biClrUsed;
		}
		else if (pbmih->biBitCount <= 8)
		{
			nOffsetBits += (ULONG)sizeof(RGBQUAD) * (1 << pbmih->biBitCount);
		}

		// Consider special cases

		if (nHeaderSize == sizeof(BITMAPINFOHEADER))
		{
			// If this is a 16 or 32 bit bitmap and BI_BITFIELDS is used, 
			// bmiColors member contains three DWORD color masks.
			// For V4 or V5 headers, this info is included in the header

			if (pbmih->biCompression == BI_BITFIELDS)
			{
				nOffsetBits += 3 * sizeof(DWORD);
			}
		}
		else if (nHeaderSize >= sizeof(BITMAPV5HEADER))
		{
			// If this is a V5 header and an ICM profile is specified,
			// we need to consider the profile data size

			PBITMAPV5HEADER pbV5h = (PBITMAPV5HEADER)pDib;

			// if the profile data comes before the pixel data, add it

			if (pbV5h->bV5ProfileData <= nOffsetBits)
			{
				nOffsetBits += pbV5h->bV5ProfileSize;
			}
		}
	}

	return nOffsetBits;
}


//////////////////////////////////////////////////////////////////////////
//
// FixBitmapHeight
//

BOOL FixBitmapHeight(PVOID pDib, ULONG nSize, BOOL bTopDown)
{
	ULONG nHeaderSize = GetBitmapHeaderSize(pDib);

	if (nHeaderSize == 0)
	{
		return FALSE;
	}

	// is this an old style BITMAPCOREHEADER?

	if (nHeaderSize == sizeof(BITMAPCOREHEADER))
	{
		PBITMAPCOREHEADER pbmch = (PBITMAPCOREHEADER)pDib;

		// fix the height value if necessary

		if (pbmch->bcHeight == 0)
		{
			// start the calculation with the header size

			ULONG nSizeImage = nSize - nHeaderSize;

			// subtract the color table size

			if (pbmch->bcBitCount <= 8)
			{
				nSizeImage -= (ULONG)sizeof(RGBTRIPLE) * (1 << pbmch->bcBitCount);
			}

			// calculate the height

			ULONG nWidth = GetBitmapLineWidthInBytes(pbmch->bcWidth, pbmch->bcBitCount);

			if (nWidth == 0)
			{
				return FALSE;
			}

			LONG nHeight = nSizeImage / nWidth;

			pbmch->bcHeight = (WORD)nHeight;
		}
	}
	else
	{
		// this is at least a BITMAPINFOHEADER

		PBITMAPINFOHEADER pbmih = (PBITMAPINFOHEADER)pDib;

		// fix the height value if necessary

		if (pbmih->biHeight == 0)
		{
			// find the size of the image data

			ULONG nSizeImage;

			if (pbmih->biSizeImage != 0)
			{
				// if the size is specified in the header, take it

				nSizeImage = pbmih->biSizeImage;
			}
			else
			{
				// start the calculation with the header size

				nSizeImage = nSize - nHeaderSize;

				// subtract the color table size

				if (pbmih->biClrUsed != 0)
				{
					nSizeImage -= sizeof(RGBQUAD) * pbmih->biClrUsed;
				}
				else if (pbmih->biBitCount <= 8)
				{
					nSizeImage -= (ULONG)sizeof(RGBQUAD) * (1 << pbmih->biBitCount);
				}

				// Consider special cases

				if (nHeaderSize == sizeof(BITMAPINFOHEADER))
				{
					// If this is a 16 or 32 bit bitmap and BI_BITFIELDS is used, 
					// bmiColors member contains three DWORD color masks.
					// For V4 or V5 headers, this info is included the header

					if (pbmih->biCompression == BI_BITFIELDS)
					{
						nSizeImage -= 3 * sizeof(DWORD);
					}
				}
				else if (nHeaderSize >= sizeof(BITMAPV5HEADER))
				{
					// If this is a V5 header and an ICM profile is specified,
					// we need to consider the profile data size

					PBITMAPV5HEADER pbV5h = (PBITMAPV5HEADER)pDib;

					// add the profile data size

					nSizeImage -= pbV5h->bV5ProfileSize;
				}

				// store the image size

				pbmih->biSizeImage = nSizeImage;
			}

			// finally, calculate the height

			ULONG nWidth = GetBitmapLineWidthInBytes(pbmih->biWidth, pbmih->biBitCount);

			if (nWidth == 0)
			{
				return FALSE;
			}

			LONG nHeight = nSizeImage / nWidth;

			pbmih->biHeight = bTopDown ? -nHeight : nHeight;
		}
	}

	return TRUE;
}


//////////////////////////////////////////////////////////////////////////
//
// FillBitmapFileHeader
//

BOOL FillBitmapFileHeader(LPCVOID pDib, PBITMAPFILEHEADER pbmfh)
{
	ULONG nSize = GetBitmapSize(pDib);

	if (nSize == 0)
	{
		return FALSE;
	}

	ULONG nOffset = GetBitmapOffsetBits(pDib);

	if (nOffset == 0)
	{
		return FALSE;
	}

	pbmfh->bfType = MAKEWORD('B', 'M');
	pbmfh->bfSize = sizeof(BITMAPFILEHEADER) + nSize;
	pbmfh->bfReserved1 = 0;
	pbmfh->bfReserved2 = 0;
	pbmfh->bfOffBits = sizeof(BITMAPFILEHEADER) + nOffset;

	return TRUE;
}

void CalcDitherPixel(int* pDitherPixel, const Color& c, const BYTE* clamp, const short* rowerr, int cursor, const bool noBias)
{
	if (noBias) {
		pDitherPixel[0] = clamp[((rowerr[cursor] + 0x1008) >> 4) + c.GetR()];
		pDitherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.GetG()];
		pDitherPixel[2] = clamp[((rowerr[cursor + 2] + 0x1008) >> 4) + c.GetB()];
		pDitherPixel[3] = clamp[((rowerr[cursor + 3] + 0x1008) >> 4) + c.GetA()];
	}
	else {
		pDitherPixel[0] = clamp[((rowerr[cursor] + 0x2010) >> 5) + c.GetR()];
		pDitherPixel[1] = clamp[((rowerr[cursor + 1] + 0x1008) >> 4) + c.GetG()];
		pDitherPixel[2] = clamp[((rowerr[cursor + 2] + 0x2010) >> 5) + c.GetB()];
		pDitherPixel[3] = c.GetA();
	}
}

bool dither_image(const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height)
{
	UINT pixelIndex = 0;

	const int DJ = 4;
	const int BLOCK_SIZE = 256;
	const int DITHER_MAX = 20;
	const int err_len = (width + 2) * DJ;
	auto clamp = make_unique <BYTE[]>(DJ * BLOCK_SIZE);
	auto erowErr = make_unique<short[]>(err_len);
	auto orowErr = make_unique<short[]>(err_len);
	auto limtb = make_unique<char[]>(2 * BLOCK_SIZE);
	auto lookup = make_unique<short[]>(65536);
	auto pDitherPixel = make_unique<int[]>(DJ);

	for (int i = 0; i < BLOCK_SIZE; ++i) {
		clamp[i] = 0;
		clamp[i + BLOCK_SIZE] = static_cast<BYTE>(i);
		clamp[i + BLOCK_SIZE * 2] = BYTE_MAX;
		clamp[i + BLOCK_SIZE * 3] = BYTE_MAX;

		limtb[i] = -DITHER_MAX;
		limtb[i + BLOCK_SIZE] = DITHER_MAX;
	}
	for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
		limtb[i + BLOCK_SIZE] = i;

	auto row0 = erowErr.get();
	auto row1 = orowErr.get();

	bool noBias = (transparentPixelIndex >= 0 || hasSemiTransparency) || nMaxColors < 64;
	int dir = 1;
	for (int i = 0; i < height; ++i) {
		if (dir < 0)
			pixelIndex += width - 1;

		int cursor0 = DJ, cursor1 = width * DJ;
		row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
		for (UINT j = 0; j < width; ++j) {
			Color c(pixels[pixelIndex]);

			CalcDitherPixel(pDitherPixel.get(), c, clamp.get(), row0, cursor0, noBias);
			int r_pix = pDitherPixel[0];
			int g_pix = pDitherPixel[1];
			int b_pix = pDitherPixel[2];
			int a_pix = pDitherPixel[3];
			auto argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
			Color c1(argb);
			if (noBias) {
				int offset = GetARGBIndex(c1, hasSemiTransparency, transparentPixelIndex >= 0);
				if (!lookup[offset])
					lookup[offset] = ditherFn(pPalette, nMaxColors, argb) + 1;
				qPixels[pixelIndex] = lookup[offset] - 1;
			}
			else
				qPixels[pixelIndex] = ditherFn(pPalette, nMaxColors, argb);

			Color c2(pPalette->Entries[qPixels[pixelIndex]]);

			r_pix = limtb[c1.GetR() - c2.GetR() + BLOCK_SIZE];
			g_pix = limtb[c1.GetG() - c2.GetG() + BLOCK_SIZE];
			b_pix = limtb[c1.GetB() - c2.GetB() + BLOCK_SIZE];
			a_pix = limtb[c1.GetA() - c2.GetA() + BLOCK_SIZE];

			int k = r_pix * 2;
			row1[cursor1 - DJ] = r_pix;
			row1[cursor1 + DJ] += (r_pix += k);
			row1[cursor1] += (r_pix += k);
			row0[cursor0 + DJ] += (r_pix + k);

			k = g_pix * 2;
			row1[cursor1 + 1 - DJ] = g_pix;
			row1[cursor1 + 1 + DJ] += (g_pix += k);
			row1[cursor1 + 1] += (g_pix += k);
			row0[cursor0 + 1 + DJ] += (g_pix + k);

			k = b_pix * 2;
			row1[cursor1 + 2 - DJ] = b_pix;
			row1[cursor1 + 2 + DJ] += (b_pix += k);
			row1[cursor1 + 2] += (b_pix += k);
			row0[cursor0 + 2 + DJ] += (b_pix + k);

			k = a_pix * 2;
			row1[cursor1 + 3 - DJ] = a_pix;
			row1[cursor1 + 3 + DJ] += (a_pix += k);
			row1[cursor1 + 3] += (a_pix += k);
			row0[cursor0 + 3 + DJ] += (a_pix + k);

			cursor0 += DJ;
			cursor1 -= DJ;
			pixelIndex += dir;
		}
		if ((i % 2) == 1)
			pixelIndex += width + 1;

		dir *= -1;
		swap(row0, row1);
	}
	return true;
}

bool dithering_image(const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, const UINT nMaxColors, ARGB* qPixels, const UINT width, const UINT height)
{
	UINT pixelIndex = 0;
	bool hasTransparency = (transparentPixelIndex >= 0 || hasSemiTransparency);
	const int DJ = 4;
	const int BLOCK_SIZE = 256;
	const int DITHER_MAX = 20;
	const int err_len = (width + 2) * DJ;
	auto clamp = make_unique <BYTE[]>(DJ * BLOCK_SIZE);
	auto erowErr = make_unique<short[]>(err_len);
	auto orowErr = make_unique<short[]>(err_len);
	auto limtb = make_unique<char[]>(2 * BLOCK_SIZE);
	auto lookup = make_unique<short[]>(65536);
	auto pDitherPixel = make_unique<int[]>(DJ);

	for (int i = 0; i < BLOCK_SIZE; ++i) {
		clamp[i] = 0;
		clamp[i + BLOCK_SIZE] = static_cast<BYTE>(i);
		clamp[i + BLOCK_SIZE * 2] = BYTE_MAX;
		clamp[i + BLOCK_SIZE * 3] = BYTE_MAX;

		limtb[i] = -DITHER_MAX;
		limtb[i + BLOCK_SIZE] = DITHER_MAX;
	}
	for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
		limtb[i + BLOCK_SIZE] = i;

	auto row0 = erowErr.get();
	auto row1 = orowErr.get();
	int dir = 1;
	for (int i = 0; i < height; ++i) {
		if (dir < 0)
			pixelIndex += width - 1;

		int cursor0 = DJ, cursor1 = width * DJ;
		row1[cursor1] = row1[cursor1 + 1] = row1[cursor1 + 2] = row1[cursor1 + 3] = 0;
		for (UINT j = 0; j < width; ++j) {
			Color c(pixels[pixelIndex]);

			CalcDitherPixel(pDitherPixel.get(), c, clamp.get(), row0, cursor0, hasTransparency);
			int r_pix = pDitherPixel[0];
			int g_pix = pDitherPixel[1];
			int b_pix = pDitherPixel[2];
			int a_pix = pDitherPixel[3];
			auto argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
			Color c1(argb);
			if (nMaxColors < 64) {
				int offset = GetARGBIndex(c1, hasSemiTransparency, transparentPixelIndex >= 0);
				if (!lookup[offset])
					lookup[offset] = ditherFn(pPalette, nMaxColors, argb) + 1;
				qPixels[pixelIndex] = lookup[offset] - 1;
			}
			else
				qPixels[pixelIndex] = ditherFn(pPalette, nMaxColors, argb);

			Color c2(pPalette->Entries[qPixels[pixelIndex]]);
			qPixels[pixelIndex] = hasSemiTransparency ? c2.GetValue() : GetARGBIndex(c2, false, transparentPixelIndex >= 0);

			r_pix = limtb[BLOCK_SIZE + c1.GetR() - c2.GetR()];
			g_pix = limtb[BLOCK_SIZE + c1.GetG() - c2.GetG()];
			b_pix = limtb[BLOCK_SIZE + c1.GetB() - c2.GetB()];
			a_pix = limtb[BLOCK_SIZE + c1.GetA() - c2.GetA()];

			int k = r_pix * 2;
			row1[cursor1 - DJ] = r_pix;
			row1[cursor1 + DJ] += (r_pix += k);
			row1[cursor1] += (r_pix += k);
			row0[cursor0 + DJ] += (r_pix + k);

			k = g_pix * 2;
			row1[cursor1 + 1 - DJ] = g_pix;
			row1[cursor1 + 1 + DJ] += (g_pix += k);
			row1[cursor1 + 1] += (g_pix += k);
			row0[cursor0 + 1 + DJ] += (g_pix + k);

			k = b_pix * 2;
			row1[cursor1 + 2 - DJ] = b_pix;
			row1[cursor1 + 2 + DJ] += (b_pix += k);
			row1[cursor1 + 2] += (b_pix += k);
			row0[cursor0 + 2 + DJ] += (b_pix + k);

			k = a_pix * 2;
			row1[cursor1 + 3 - DJ] = a_pix;
			row1[cursor1 + 3 + DJ] += (a_pix += k);
			row1[cursor1 + 3] += (a_pix += k);
			row0[cursor0 + 3 + DJ] += (a_pix + k);

			cursor0 += DJ;
			cursor1 -= DJ;
			pixelIndex += dir;
		}
		if ((i % 2) == 1)
			pixelIndex += (width + 1);

		dir *= -1;
		swap(row0, row1);
	}
	return true;
}

bool ProcessImagePixels(Bitmap* pDest, const ARGB* qPixels, const bool& hasSemiTransparency, const int& transparentPixelIndex)
{
	UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
	if (bpp < 16)
		return false;

	if (hasSemiTransparency && pDest->GetPixelFormat() < PixelFormat32bppARGB)
		pDest->ConvertFormat(PixelFormat32bppARGB, DitherTypeNone, PaletteTypeCustom, nullptr, 0);
	else if (transparentPixelIndex >= 0 && pDest->GetPixelFormat() < PixelFormat16bppARGB1555)
		pDest->ConvertFormat(PixelFormat16bppARGB1555, DitherTypeNone, PaletteTypeCustom, nullptr, 0);
	else if (pDest->GetPixelFormat() != PixelFormat16bppRGB565)
		pDest->ConvertFormat(PixelFormat16bppRGB565, DitherTypeNone, PaletteTypeCustom, nullptr, 0);

	BitmapData targetData;
	UINT w = pDest->GetWidth();
	UINT h = pDest->GetHeight();

	Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
	if (status != Ok) {
		cerr << "Cannot write image" << endl;
		return false;
	}

	int pixelIndex = 0;

	auto pRowDest = (LPBYTE)targetData.Scan0;
	UINT strideDest;

	// Compensate for possible negative stride
	if (targetData.Stride > 0)
		strideDest = targetData.Stride;
	else {
		pRowDest += h * targetData.Stride;
		strideDest = -targetData.Stride;
	}

	bpp = GetPixelFormatSize(pDest->GetPixelFormat());
	if (bpp == 32) {
		for (UINT y = 0; y < h; ++y) {	// For each row...
			for (UINT x = 0; x < w * 4;) {
				Color c(qPixels[pixelIndex++]);
				pRowDest[x++] = c.GetB();
				pRowDest[x++] = c.GetG();
				pRowDest[x++] = c.GetR();
				pRowDest[x++] = c.GetA();
			}
			pRowDest += strideDest;
		}
	}
	else if (bpp == 16) {
		for (UINT y = 0; y < h; ++y) {	// For each row...
			for (UINT x = 0; x < w * 2;) {
				auto argb = qPixels[pixelIndex++];
				pRowDest[x++] = static_cast<BYTE>(argb & 0xFF);
				pRowDest[x++] = static_cast<BYTE>(argb >> 8);
			}
			pRowDest += strideDest;
		}
	}
	else {
		for (UINT y = 0; y < h; ++y) {	// For each row...
			for (UINT x = 0; x < w; ++x) {
				auto argb = qPixels[pixelIndex++];
				pRowDest[x] = static_cast<BYTE>(argb);
			}
			pRowDest += strideDest;
		}
	}

	status = pDest->UnlockBits(&targetData);
	return pDest->GetLastStatus() == Ok;
}

bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, const unsigned short* qPixels, const bool hasTransparent)
{
	if (hasTransparent) {
		BYTE value = 0;

		auto pPropertyItem = make_unique<PropertyItem>();
		pPropertyItem.get()->id = PropertyTagIndexTransparent;
		pPropertyItem.get()->length = 1;
		pPropertyItem.get()->type = PropertyTagTypeByte;
		pPropertyItem.get()->value = &value;

		pDest->SetPropertyItem(pPropertyItem.get());
	}

	pDest->SetPalette(pPalette);

	BitmapData targetData;
	UINT w = pDest->GetWidth();
	UINT h = pDest->GetHeight();

	Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
	if (status != Ok) {
		cerr << "Cannot write image" << endl;
		return false;
	}

	int pixelIndex = 0;

	auto pRowDest = (LPBYTE)targetData.Scan0;
	UINT strideDest;

	// Compensate for possible negative stride
	if (targetData.Stride > 0)
		strideDest = targetData.Stride;
	else {
		pRowDest += h * targetData.Stride;
		strideDest = -targetData.Stride;
	}

	UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
	// Second loop: fill indexed bitmap
	for (UINT y = 0; y < h; ++y) {	// For each row...
		for (UINT x = 0; x < w; ++x) {	// ...for each pixel...
			BYTE nibbles = 0;
			BYTE index = static_cast<BYTE>(qPixels[pixelIndex++]);

			switch (bpp)
			{
			case 8:
				pRowDest[x] = index;
				break;
			case 4:
				// First pixel is the high nibble. From and To indices are 0..16
				nibbles = pRowDest[x / 2];
				if ((x & 1) == 0) {
					nibbles &= 0x0F;
					nibbles |= (BYTE)(index << 4);
				}
				else {
					nibbles &= 0xF0;
					nibbles |= index;
				}

				pRowDest[x / 2] = nibbles;
				break;
			case 1:
				// First pixel is MSB. From and To are 0 or 1.
				int pos = x / 8;
				BYTE mask = (BYTE)(128 >> (x & 7));
				if (index == 0)
					pRowDest[pos] &= (BYTE)~mask;
				else
					pRowDest[pos] |= mask;
				break;
			}
		}

		pRowDest += strideDest;
	}

	status = pDest->UnlockBits(&targetData);
	return pDest->GetLastStatus() == Ok;
}

bool GrabPixels(Bitmap* pSource, vector<ARGB>& pixels, bool& hasSemiTransparency, int& transparentPixelIndex, ARGB& transparentColor)
{
	const UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
	const UINT bitmapWidth = pSource->GetWidth();
	const UINT bitmapHeight = pSource->GetHeight();

	hasSemiTransparency = false;
	transparentPixelIndex = -1;

	int transparentIndex = -1;
	int paletteSize = pSource->GetPaletteSize();
	auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + (1 << bitDepth) * sizeof(ARGB));
	auto pPalette = (ColorPalette*)pPaletteBytes.get();
	if (paletteSize > 0)
		pSource->GetPalette(pPalette, paletteSize);

	auto nSize = pSource->GetPropertyItemSize(PropertyTagIndexTransparent);
	if (nSize > 0) {
		auto pPropertyItem = make_unique<PropertyItem[]>(nSize);
		pSource->GetPropertyItem(PropertyTagIndexTransparent, nSize, pPropertyItem.get());
		if (pPropertyItem.get()->length > 0) {
			transparentIndex = *(BYTE*)pPropertyItem.get()->value;
			Color c(pPalette->Entries[transparentIndex]);
			pPalette->Entries[transparentIndex] = transparentColor = Color::MakeARGB(0, c.GetR(), c.GetG(), c.GetB());
		}
	}

	int pixelIndex = 0;
	if (pSource->GetPixelFormat() == PixelFormat8bppIndexed) {
		BitmapData data;
		Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, pSource->GetPixelFormat(), &data);
		if (status != Ok)
			return false;

		auto pRowSource = (LPBYTE)data.Scan0;
		UINT strideSource;

		if (data.Stride > 0)
			strideSource = data.Stride;
		else {
			// Compensate for possible negative stride
			// (not needed for first loop, but we have to do it
			// for second loop anyway)
			pRowSource += bitmapHeight * data.Stride;
			strideSource = -data.Stride;
		}

		// First loop: gather color information
		for (UINT y = 0; y < bitmapHeight; ++y) {	// For each row...
			auto pPixelSource = pRowSource;

			for (UINT x = 0; x < bitmapWidth; ++x) {	// ...for each pixel...
				BYTE index = *pPixelSource++;
				auto argb = pPalette->Entries[index];
				Color c(argb);

				if (c.GetA() < BYTE_MAX) {
					if (c.GetA() == 0) {
						transparentColor = argb;
						transparentPixelIndex = pixelIndex;
						if (transparentColor == 0 && transparentIndex < 0)
							argb = transparentColor = Color::MakeARGB(0, 51, 102, 102);
					}
					else
						hasSemiTransparency = true;
				}
				pixels[pixelIndex++] = argb;
			}

			pRowSource += strideSource;
		}

		pSource->UnlockBits(&data);

		return true;
	}

	BitmapData data;
	Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, PixelFormat32bppARGB, &data);
	if (status != Ok)
		return false;

	auto pRowSource = (LPBYTE)data.Scan0;
	UINT strideSource;

	if (data.Stride > 0) strideSource = data.Stride;
	else
	{
		// Compensate for possible negative stride
		// (not needed for first loop, but we have to do it
		// for second loop anyway)
		pRowSource += bitmapHeight * data.Stride;
		strideSource = -data.Stride;
	}

	// First loop: gather color information
	for (UINT y = 0; y < bitmapHeight; ++y) {	// For each row...
		auto pPixelSource = pRowSource;

		for (UINT x = 0; x < bitmapWidth; ++x) {	// ...for each pixel...
			BYTE pixelBlue = *pPixelSource++;
			BYTE pixelGreen = *pPixelSource++;
			BYTE pixelRed = *pPixelSource++;
			BYTE pixelAlpha = *pPixelSource++;
			auto argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
			if (transparentIndex > -1 && Color(transparentColor).ToCOLORREF() == Color(argb).ToCOLORREF()) {
				pixelAlpha = 0;
				argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
			}

			if (pixelAlpha < BYTE_MAX) {
				if (pixelAlpha == 0) {
					transparentColor = argb;
					transparentPixelIndex = pixelIndex;
					if (transparentColor == 0 && transparentIndex < 0)
						argb = transparentColor = Color::MakeARGB(0, 51, 102, 102);
				}
				else
					hasSemiTransparency = true;
			}
			pixels[pixelIndex++] = argb;
		}

		pRowSource += strideSource;
	}

	pSource->UnlockBits(&data);

	return true;
}

bool HasTransparency(Bitmap* pSource)
{
	const UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
	const UINT bitmapWidth = pSource->GetWidth();
	const UINT bitmapHeight = pSource->GetHeight();

	// Not an alpha-capable color format. Note that GDI+ indexed images are alpha-capable on the palette.
	if (!(pSource->GetFlags() & ImageFlags::ImageFlagsHasAlpha))
		return false;

	if (pSource->GetPixelFormat() & PixelFormatIndexed) {
		int paletteSize = pSource->GetPaletteSize();
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + (1 << bitDepth) * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pSource->GetPalette(pPalette, paletteSize);

		for (UINT i = 0; i < pPalette->Count; ++i) {
			Color c(pPalette->Entries[i]);

			if (c.GetA() < BYTE_MAX)
				return true;
		}

		auto nSize = pSource->GetPropertyItemSize(PropertyTagIndexTransparent);
		if (nSize > 0) {
			auto pPropertyItem = make_unique<PropertyItem[]>(nSize);
			pSource->GetPropertyItem(PropertyTagIndexTransparent, nSize, pPropertyItem.get());
			if (pPropertyItem.get()->length > 0)
				return true;
		}
		return false;
	}

	BitmapData data;
	Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, PixelFormat32bppARGB, &data);
	if (status != Ok)
		return false;

	auto pRowSource = (LPBYTE)data.Scan0;
	UINT strideSource;

	if (data.Stride > 0) strideSource = data.Stride;
	else
	{
		// Compensate for possible negative stride
		// (not needed for first loop, but we have to do it
		// for second loop anyway)
		pRowSource += bitmapHeight * data.Stride;
		strideSource = -data.Stride;
	}

	// First loop: gather color information
	for (UINT y = 0; y < bitmapHeight; ++y) {	// For each row...
		auto pPixelSource = pRowSource;

		for (UINT x = 0; x < bitmapWidth; ++x) {	// ...for each pixel...
			BYTE pixelBlue = *pPixelSource++;
			BYTE pixelGreen = *pPixelSource++;
			BYTE pixelRed = *pPixelSource++;
			BYTE pixelAlpha = *pPixelSource++;

			if (pixelAlpha < BYTE_MAX) {
				pSource->UnlockBits(&data);
				return true;
			}
		}

		pRowSource += strideSource;
	}

	pSource->UnlockBits(&data);

	return false;
}