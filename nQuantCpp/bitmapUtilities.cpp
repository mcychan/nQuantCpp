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

BOOL GetBitmapDimensions(LPCVOID pDib, UINT *pWidth, UINT *pHeight)
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

bool dither_image(const ARGB* pixels, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, short* qPixels, const UINT width, const UINT height)
{
	UINT pixelIndex = 0;
	bool odd_scanline = false;
	short *row0, *row1;
	int a_pix, r_pix, g_pix, b_pix, dir, k;
	const int DJ = 4;
	const int DITHER_MAX = 20;
	const int err_len = (width + 2) * DJ;
	byte clamp[DJ * 256] = { 0 };
	auto erowErr = make_unique<short[]>(err_len);
	auto orowErr = make_unique<short[]>(err_len);
	char limtb[512] = { 0 };
	auto lim = &limtb[256];
	auto erowerr = erowErr.get();
	auto orowerr = orowErr.get();
	int lookup[65536] = { 0 };

	for (int i = 0; i < 256; i++) {
		clamp[i] = 0;
		clamp[i + 256] = static_cast<byte>(i);
		clamp[i + 512] = BYTE_MAX;
		clamp[i + 768] = BYTE_MAX;

		limtb[i] = -DITHER_MAX;
		limtb[i + 256] = DITHER_MAX;
	}
	for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
		limtb[i + 256] = i;

	for (int i = 0; i < height; i++) {
		if (odd_scanline) {
			dir = -1;
			pixelIndex += (width - 1);
			row0 = &orowerr[DJ];
			row1 = &erowerr[width * DJ];
		}
		else {
			dir = 1;
			row0 = &erowerr[DJ];
			row1 = &orowerr[width * DJ];
		}
		row1[0] = row1[1] = row1[2] = row1[3] = 0;
		for (UINT j = 0; j < width; j++) {
			Color c(pixels[pixelIndex]);

			r_pix = clamp[((row0[0] + 0x1008) >> 4) + c.GetR()];
			g_pix = clamp[((row0[1] + 0x1008) >> 4) + c.GetG()];
			b_pix = clamp[((row0[2] + 0x1008) >> 4) + c.GetB()];
			a_pix = clamp[((row0[3] + 0x1008) >> 4) + c.GetA()];

			ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
			Color c1(argb);
			int offset = getARGBIndex(c1, hasSemiTransparency, transparentPixelIndex);
			if (!lookup[offset]) {
				auto argb1 = Color::MakeARGB(BYTE_MAX, (c1.GetR() & 0xF8), (c1.GetG() & 0xFC), (c1.GetB() & 0xF8));
				if (hasSemiTransparency)
					argb1 = Color::MakeARGB((c1.GetA() & 0xF0), (c1.GetR() & 0xF0), (c1.GetG() & 0xF0), (c1.GetB() & 0xF0));
				else if (transparentPixelIndex >= 0)
					argb1 = Color::MakeARGB((c1.GetA() < BYTE_MAX) ? 0 : BYTE_MAX, (c1.GetR() & 0xF8), (c1.GetG() & 0xF8), (c1.GetB() & 0xF8));
				lookup[offset] = argb1 + 1;
			}
			qPixels[pixelIndex] = static_cast<short>(offset);

			Color c2(static_cast<ARGB>(lookup[offset] - 1));

			r_pix = lim[r_pix - c2.GetR()];
			g_pix = lim[g_pix - c2.GetG()];
			b_pix = lim[b_pix - c2.GetB()];
			a_pix = lim[a_pix - c2.GetA()];

			k = r_pix * 2;
			row1[0 - DJ] = r_pix;
			row1[0 + DJ] += (r_pix += k);
			row1[0] += (r_pix += k);
			row0[0 + DJ] += (r_pix += k);

			k = g_pix * 2;
			row1[1 - DJ] = g_pix;
			row1[1 + DJ] += (g_pix += k);
			row1[1] += (g_pix += k);
			row0[1 + DJ] += (g_pix += k);

			k = b_pix * 2;
			row1[2 - DJ] = b_pix;
			row1[2 + DJ] += (b_pix += k);
			row1[2] += (b_pix += k);
			row0[2 + DJ] += (b_pix += k);

			k = a_pix * 2;
			row1[3 - DJ] = a_pix;
			row1[3 + DJ] += (a_pix += k);
			row1[3] += (a_pix += k);
			row0[3 + DJ] += (a_pix += k);

			row0 += DJ;
			row1 -= DJ;
			pixelIndex += dir;
		}
		if ((i % 2) == 1)
			pixelIndex += (width + 1);

		odd_scanline = !odd_scanline;
	}
	return true;
}

bool dither_image(const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, const UINT nMaxColors, short* qPixels, const UINT width, const UINT height)
{
	UINT pixelIndex = 0;
	
	bool odd_scanline = false;
	short *row0, *row1;
	int j, a_pix, r_pix, g_pix, b_pix, dir, k;
	const int DJ = 4;
	const int DITHER_MAX = 20;
	const int err_len = (width + 2) * DJ;
	byte clamp[DJ * 256] = { 0 };
	auto erowErr = make_unique<short[]>(err_len);
	auto orowErr = make_unique<short[]>(err_len);
	char limtb[512] = { 0 };
	auto lim = &limtb[256];
	auto erowerr = erowErr.get();
	auto orowerr = orowErr.get();
	short lookup[65536] = { 0 };

	for (int i = 0; i < 256; i++) {
		clamp[i] = 0;
		clamp[i + 256] = static_cast<byte>(i);
		clamp[i + 512] = BYTE_MAX;
		clamp[i + 768] = BYTE_MAX;

		limtb[i] = -DITHER_MAX;
		limtb[i + 256] = DITHER_MAX;
	}
	for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
		limtb[i + 256] = i;

	for (int i = 0; i < height; i++) {
		if (odd_scanline) {
			dir = -1;
			pixelIndex += (width - 1);
			row0 = &orowerr[DJ];
			row1 = &erowerr[width * DJ];
		}
		else {
			dir = 1;
			row0 = &erowerr[DJ];
			row1 = &orowerr[width * DJ];
		}
		row1[0] = row1[1] = row1[2] = row1[3] = 0;
		for (j = 0; j < width; j++) {
			Color c(pixels[pixelIndex]);

			r_pix = clamp[((row0[0] + 0x1008) >> 4) + c.GetR()];
			g_pix = clamp[((row0[1] + 0x1008) >> 4) + c.GetG()];
			b_pix = clamp[((row0[2] + 0x1008) >> 4) + c.GetB()];
			a_pix = clamp[((row0[3] + 0x1008) >> 4) + c.GetA()];

			ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
			Color c1(argb);
			int offset = getARGBIndex(c1, hasSemiTransparency, transparentPixelIndex);
			if (!lookup[offset])
				lookup[offset] = ditherFn(pPalette, nMaxColors, c.GetA() ? argb : pixels[pixelIndex]) + 1;
			qPixels[pixelIndex] = lookup[offset] - 1;

			Color c2(pPalette->Entries[qPixels[pixelIndex]]);

			r_pix = lim[r_pix - c2.GetR()];
			g_pix = lim[g_pix - c2.GetG()];
			b_pix = lim[b_pix - c2.GetB()];
			a_pix = lim[a_pix - c2.GetA()];

			k = r_pix * 2;
			row1[0 - DJ] = r_pix;
			row1[0 + DJ] += (r_pix += k);
			row1[0] += (r_pix += k);
			row0[0 + DJ] += (r_pix += k);

			k = g_pix * 2;
			row1[1 - DJ] = g_pix;
			row1[1 + DJ] += (g_pix += k);
			row1[1] += (g_pix += k);
			row0[1 + DJ] += (g_pix += k);

			k = b_pix * 2;
			row1[2 - DJ] = b_pix;
			row1[2 + DJ] += (b_pix += k);
			row1[2] += (b_pix += k);
			row0[2 + DJ] += (b_pix += k);

			k = a_pix * 2;
			row1[3 - DJ] = a_pix;
			row1[3 + DJ] += (a_pix += k);
			row1[3] += (a_pix += k);
			row0[3 + DJ] += (a_pix += k);

			row0 += DJ;
			row1 -= DJ;
			pixelIndex += dir;
		}
		if ((i % 2) == 1)
			pixelIndex += (width + 1);

		odd_scanline = !odd_scanline;
	}
	return true;
}

bool ProcessImagePixels(Bitmap* pDest, const short* qPixels)
{
	UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
	if (bpp < 16)
		return false;

	BitmapData targetData;
	UINT w = pDest->GetWidth();
	UINT h = pDest->GetHeight();

	Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
	if (status != Ok) {
		cerr << "Cannot write image" << endl;
		return false;
	}

	int pixelIndex = 0;

	auto pRowDest = (byte*)targetData.Scan0;
	UINT strideDest;

	// Compensate for possible negative stride
	if (targetData.Stride > 0)
		strideDest = targetData.Stride;
	else {
		pRowDest += h * targetData.Stride;
		strideDest = -targetData.Stride;
	}

	// Second loop: fill indexed bitmap
	for (UINT y = 0; y < h; y++) {	// For each row...
		for (UINT x = 0; x < w * 2;) {
			auto argb = static_cast<unsigned short>(qPixels[pixelIndex++]);
			pRowDest[x++] = static_cast<byte>(argb & 0xFF);
			pRowDest[x++] = static_cast<byte>(argb >> 8);
		}
		pRowDest += strideDest;
	}

	status = pDest->UnlockBits(&targetData);
	return pDest->GetLastStatus() == Ok;
}

bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, const short* qPixels)
{
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

	auto pRowDest = (byte*)targetData.Scan0;
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
	for (UINT y = 0; y < h; y++) {	// For each row...
		for (UINT x = 0; x < w; x++) {	// ...for each pixel...
			byte nibbles = 0;
			byte index = static_cast<byte>(qPixels[pixelIndex++]);

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
					nibbles |= (byte)(index << 4);
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
				byte mask = (byte)(128 >> (x & 7));
				if (index == 0)
					pRowDest[pos] &= (byte)~mask;
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
		
	int pixelIndex = 0;
	if (bitDepth <= 16) {
		for (UINT y = 0; y < bitmapHeight; y++) {
			for (UINT x = 0; x < bitmapWidth; x++) {
				Color color;
				pSource->GetPixel(x, y, &color);
				if (color.GetA() < BYTE_MAX) {
					hasSemiTransparency = true;
					if (color.GetA() == 0) {
						transparentPixelIndex = pixelIndex;
						transparentColor = color.GetValue();
					}
				}
				pixels[pixelIndex++] = color.GetValue();
			}
		}
	}

	// Lock bits on 3x8 source bitmap
	else {
		BitmapData data;
		Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, pSource->GetPixelFormat(), &data);
		if (status != Ok)
			return false;

		auto pRowSource = (byte*)data.Scan0;
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
		for (UINT y = 0; y < bitmapHeight; y++) {	// For each row...
			auto pPixelSource = pRowSource;

			for (UINT x = 0; x < bitmapWidth; x++) {	// ...for each pixel...
				byte pixelBlue = *pPixelSource++;
				byte pixelGreen = *pPixelSource++;
				byte pixelRed = *pPixelSource++;
				byte pixelAlpha = bitDepth < 32 ? BYTE_MAX : *pPixelSource++;

				auto argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
				if (pixelAlpha < BYTE_MAX) {
					transparentPixelIndex = pixelIndex;
					if (pixelAlpha == 0)
						transparentColor = argb;
				}
				pixels[pixelIndex++] = argb;
			}

			pRowSource += strideSource;
		}

		pSource->UnlockBits(&data);
	}
	return true;
}
