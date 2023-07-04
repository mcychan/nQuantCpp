#pragma once
#include "CIELABConvertor.h"
#include <memory>
#include <unordered_map>
#include <vector>
using namespace std;

namespace PnnLABQuant
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

	class PnnLABQuantizer
	{
		private:
			bool hasSemiTransparency = false;
			int m_transparentPixelIndex = -1;
			bool isGA = false;
			double proportional = 1.0, ratio = .5, ratioY = .5;
			unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
			unordered_map<ARGB, vector<unsigned short> > closestMap;
			unordered_map<ARGB, unsigned short> nearestMap;
			vector<float> saliencies;

			struct pnnbin {
				float ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
				float cnt = 0;
				int nn = 0, fw = 0, bk = 0, tm = 0, mtm = 0;
			};

			void find_nn(pnnbin* bins, int idx, bool texicab);
			unsigned short closestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos);
			bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither);
			void clear();

		public:
			PnnLABQuantizer();
			PnnLABQuantizer(const PnnLABQuantizer& quantizer);
			void pnnquan(const vector<ARGB>& pixels, ColorPalette* pPalette, UINT& nMaxColors);
			bool IsGA() const;
			void getLab(const Color& c, CIELABConvertor::Lab& lab1);
			bool hasAlpha() const;
			unsigned short nearestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos);
			void setRatio(double ratioX, double ratioY);
			void grabPixels(Bitmap* srcImg, vector<ARGB>& pixels, UINT& nMaxColors, bool& hasSemiTransparency);
			bool QuantizeImage(const vector<ARGB>& pixels, const UINT bitmapWidth, ColorPalette* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither = true);
			bool QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither = true);
	};
}