#pragma once
#include "bitmapUtilities.h"

namespace Peano
{
	class GilbertCurve
	{
		public:
			static void dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, float* saliencies, double weight = 1.0);
	};
}
