#pragma once
#include "bitmapUtilities.h"

namespace BlueNoise
{
	void dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels);
}
