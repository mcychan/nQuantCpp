#pragma once
#include "bitmapUtilities.h"

namespace Peano
{
	class GilbertCurve
	{
		public:
			static void dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColor, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, float* saliencies, double weight = 1.0);

			static void dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColor, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, ARGB* qPixels, float* saliencies, double weight = 1.0);
	};
}
