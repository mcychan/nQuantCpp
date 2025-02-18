#pragma once
#include <functional>
#include <iostream>
#include <memory>
#include <vector>
using namespace std;

#ifndef _WIN32
	#include <guiddef.h>
	#include <gdiplus.h>
	using namespace Gdiplus;
#endif


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapHeaderSize
//

ULONG GetBitmapHeaderSize(LPCVOID pDib);


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapLineWidthInBytes
//

ULONG GetBitmapLineWidthInBytes(ULONG nWidthInPixels, ULONG nBitCount);

//////////////////////////////////////////////////////////////////////////
//
// GetBitmapDimensions
//

BOOL GetBitmapDimensions(LPCVOID pDib, UINT *pWidth, UINT *pHeight);


//////////////////////////////////////////////////////////////////////////
//
// GetBitmapSize
//

ULONG GetBitmapSize(LPCVOID pDib);



//////////////////////////////////////////////////////////////////////////
//
// GetBitmapOffsetBits
//

ULONG GetBitmapOffsetBits(LPCVOID pDib);

//////////////////////////////////////////////////////////////////////////
//
// FixBitmapHeight
//

BOOL FixBitmapHeight(PVOID pDib, ULONG nSize, BOOL bTopDown);


//////////////////////////////////////////////////////////////////////////
//
// FillBitmapFileHeader
//

BOOL FillBitmapFileHeader(LPCVOID pDib, PBITMAPFILEHEADER pbmfh);

using DitherFn = function<unsigned short(const ARGB*, const UINT nMaxColor, ARGB, const UINT)>;

using GetColorIndexFn = function<int(const Color&)>;

void CalcDitherPixel(int* pDitherPixel, const Color& c, const BYTE* clamp, const short* rowerr, int cursor, const bool noBias);

bool dither_image(const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColors, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, unsigned short* qPixels, const UINT width, const UINT height);

bool dithering_image(const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, const bool& hasSemiTransparency, const int& transparentPixelIndex, const UINT nMaxColors, ARGB* qPixels, const UINT width, const UINT height);

bool ProcessImagePixels(Bitmap* pDest, const ARGB* qPixels, const bool& hasSemiTransparency, const int& transparentPixelIndex);

bool ProcessImagePixels(Bitmap* pDest, const unsigned short* qPixels, const bool hasTransparent);

bool GrabPixels(Bitmap* pSource, vector<ARGB>& pixels, int& semiTransCount, int& transparentPixelIndex, ARGB& transparentColor, const BYTE alphaThreshold, const UINT nMaxColors = 2);

int GrabPixels(Bitmap* pSource, vector<ARGB>& pixels, bool& hasSemiTransparency, int& transparentPixelIndex, ARGB& transparentColor, const BYTE alphaThreshold, const UINT nMaxColors = 2);

bool HasTransparency(Bitmap* pSource);

inline int GetARGB1555(const Color& c)
{
	return (c.GetA() & 0x80) << 8 | (c.GetR() & 0xF8) << 7 | (c.GetG() & 0xF8) << 2 | (c.GetB() >> 3);
}

inline int GetARGBIndex(const Color& c, const bool hasSemiTransparency, const bool hasTransparency)
{
	if (hasSemiTransparency)
		return (c.GetA() & 0xF0) << 8 | (c.GetR() & 0xF0) << 4 | (c.GetG() & 0xF0) | (c.GetB() >> 4);
	if (hasTransparency)
		return GetARGB1555(c);
	return (c.GetR() & 0xF8) << 8 | (c.GetG() & 0xFC) << 3 | (c.GetB() >> 3);
}
