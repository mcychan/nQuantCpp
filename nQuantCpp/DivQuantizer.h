#pragma once
#include "bitmapUtilities.h"

#include <memory>
#include <vector>
using namespace std;

namespace DivQuant
{
	// =============================================================
	// Quantizer objects and functions
	//
	// COVERED CODE IS PROVIDED UNDER THIS LICENSE ON AN "AS IS" BASIS, WITHOUT WARRANTY
	// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
	// THAT THE COVERED CODE IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE
	// OR NON-INFRINGING. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE COVERED
	// CODE IS WITH YOU. SHOULD ANY COVERED CODE PROVE DEFECTIVE IN ANY RESPECT, YOU (NOT
	// THE INITIAL DEVELOPER OR ANY OTHER CONTRIBUTOR) ASSUME THE COST OF ANY NECESSARY
	// SERVICING, REPAIR OR CORRECTION. THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL
	// PART OF THIS LICENSE. NO USE OF ANY COVERED CODE IS AUTHORIZED HEREUNDER EXCEPT UNDER
	// THIS DISCLAIMER.
	//
	// Use at your own risk!
	// =============================================================

	class DivQuantizer
	{
		public:
			void quant_varpart_fast(const ARGB* inPixels, const UINT numPixels, ColorPalette* pPalette,
				const UINT numRows = 1, const bool allPixelsUnique = true,
				const int num_bits = 8, const int dec_factor = 1, const int max_iters = 10);
			bool QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither = true);
	};
}