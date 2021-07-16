#pragma once
#include "bitmapUtilities.h"

namespace Riemersma
{
	class HilbertCurve
	{
		public:
			static void dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, unsigned short* qPixels);
	};
}
