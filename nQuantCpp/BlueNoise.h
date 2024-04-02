#pragma once
#include "bitmapUtilities.h"

namespace BlueNoise
{
	extern const char TELL_BLUE_NOISE[];
	
	ARGB diffuse(const Color& pixel, const Color& qPixel, const float weight, const float strength, const int x, const int y);

	void dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const unsigned short nMaxColors, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, const float weight = 1.0f);
}
