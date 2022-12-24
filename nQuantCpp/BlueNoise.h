#pragma once
#include "bitmapUtilities.h"

namespace BlueNoise
{
	extern const char RAW_BLUE_NOISE[];
	
	ARGB diffuse(const Color& pixel, const Color& qPixel, const float weight, const float strength, const int x, const int y);

	void dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, const float weight = 1.0f);
}
