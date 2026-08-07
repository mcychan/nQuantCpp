#pragma once
#include "bitmapUtilities.h"

namespace BlueNoise
{
	extern const char TELL_BLUE_NOISE[];
	
	ARGB diffuse(const Color& pixel, const Color& qPixel, const float weight, const float strength, const int x, const int y);

	void dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const unsigned short nMaxColors, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, const float weight = 1.0f);

	ARGB dither_pixel(const Color& pixel, const float saliency, const int x, const int y,
		const float baseSpread, unsigned int frameIndex);

	bool dither_image(const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColors, DitherFn ditherFn,
		unsigned short* qPixels, const UINT width, const UINT height,
		const vector<float>& saliencies, unsigned int frameIndex);

	bool dither_image(const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColors, DitherFn ditherFn,
		unsigned short* qPixels, const UINT width, const UINT height, unsigned int frameIndex);
}
