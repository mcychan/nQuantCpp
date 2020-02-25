#pragma once
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2019 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnQuantizer.h"
#include "bitmapUtilities.h"
#include <unordered_map>

namespace PnnQuant
{
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, vector<unsigned short> > closestMap;

	struct pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		int cnt = 0;
		int nn = 0, fw = 0, bk = 0, tm = 0, mtm = 0;
	};

	void find_nn(pnnbin* bins, int idx)
	{
		int i, nn = 0;
		double err = 1e100;

		auto& bin1 = bins[idx];
		auto n1 = bin1.cnt;
		auto wa = bin1.ac;
		auto wr = bin1.rc;
		auto wg = bin1.gc;
		auto wb = bin1.bc;
		for (i = bin1.fw; i; i = bins[i].fw) {
			double nerr = sqr(bins[i].rc - wr) + sqr(bins[i].gc - wg) + sqr(bins[i].bc - wb);
			if(hasSemiTransparency)
				nerr += sqr(bins[i].ac - wa);
			double n2 = bins[i].cnt;
			nerr *= (n1 * n2) / (n1 + n2);
			if (nerr >= err)
				continue;
			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	int pnnquan(const vector<ARGB>& pixels, ColorPalette* pPalette, UINT nMaxColors, bool quan_sqrt)
	{
		auto bins = make_unique<pnnbin[]>(65536);
		auto heap = make_unique<int[]>(65537);
		double err, n1, n2;

		/* Build histogram */
		for (const auto& pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			Color c(pixel);
			int index = getARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex);
			auto& tb = bins[index];
			if (hasSemiTransparency)
				tb.ac += c.GetA();
			tb.rc += c.GetR();
			tb.gc += c.GetG();
			tb.bc += c.GetB();
			tb.cnt++;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < 65536; ++i) {
			if (!bins[i].cnt)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			if (hasSemiTransparency)
				bins[i].ac *= d;
			bins[i].rc *= d;
			bins[i].gc *= d;
			bins[i].bc *= d;
			if(quan_sqrt)
				bins[i].cnt = _sqrt(bins[i].cnt);
			bins[maxbins] = bins[i];
			++maxbins;
		}

		for (int i = 0; i < maxbins - 1; ++i) {
			bins[i].fw = i + 1;
			bins[i + 1].bk = i;
		}

		//	bins[0].bk = bins[i].fw = 0;

		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins.get(), i);
			/* Push slot on heap */
			err = bins[i].err;
			for (l = ++heap[0]; l > 1; l = l2) {
				l2 = l >> 1;
				if (bins[h = heap[l2]].err <= err)
					break;
				heap[l] = h;
			}
			heap[l] = i;
		}

		/* Merge bins which increase error the least */
		int extbins = maxbins - nMaxColors;
		for (int i = 0; i < extbins; ) {
			int b1;
			
			/* Use heap to find which bins to merge */
			for (;;) {
				auto& tb = bins[b1 = heap[1]]; /* One with least error */
											   /* Is stored error up to date? */
				if ((tb.tm >= tb.mtm) && (bins[tb.nn].mtm <= tb.tm))
					break;
				if (tb.mtm == 0xFFFF) /* Deleted node */
					b1 = heap[1] = heap[heap[0]--];
				else /* Too old error value */
				{
					find_nn(bins.get(), b1);
					tb.tm = i;
				}
				/* Push slot down */
				err = bins[b1].err;
				for (l = 1; (l2 = l + l) <= heap[0]; l = l2) {
					if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
						l2++;
					if (err <= bins[h = heap[l2]].err)
						break;
					heap[l] = h;
				}
				heap[l] = b1;
			}

			/* Do a merge */
			auto& tb = bins[b1];
			auto& nb = bins[tb.nn];
			n1 = tb.cnt;
			n2 = nb.cnt;
			double d = 1.0 / (n1 + n2);
			if (hasSemiTransparency)
				tb.ac = d * rint(n1 * tb.ac + n2 * nb.ac);
			tb.rc = d * rint(n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * rint(n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * rint(n1 * tb.bc + n2 * nb.bc);
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		UINT k = 0;
		for (int i = 0;; ++k) {
			auto alpha = hasSemiTransparency ? rint(bins[i].ac) : BYTE_MAX;
			pPalette->Entries[k] = Color::MakeARGB(alpha, rint(bins[i].rc), rint(bins[i].gc), rint(bins[i].bc));
			if (m_transparentPixelIndex >= 0 && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);

			if (!(i = bins[i].fw))
				break;
		}

		return 0;
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		unsigned short k = 0;
		Color c(argb);

		double mindist = INT_MAX;
		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette->Entries[i]);
			double curdist = sqr(c2.GetA() - c.GetA());
			if (curdist > mindist)
				continue;

			curdist += sqr(c2.GetR() - c.GetR());
			if (curdist > mindist)
				continue;

			curdist += sqr(c2.GetG() - c.GetG());
			if (curdist > mindist)
				continue;

			curdist += sqr(c2.GetB() - c.GetB());
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		return k;
	}

	unsigned short closestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		UINT k = 0;
		Color c(argb);
		vector<unsigned short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			for (; k < nMaxColors; ++k) {
				Color c2(pPalette->Entries[k]);
				closest[4] = abs(c.GetA() - c2.GetA()) + abs(c.GetR() - c2.GetR()) + abs(c.GetG() - c2.GetG()) + abs(c.GetB() - c2.GetB());
				if (closest[4] < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = closest[4];
				}
				else if (closest[4] < closest[3]) {
					closest[1] = k;
					closest[3] = closest[4];
				}
			}

			if (closest[3] == SHORT_MAX)
				closest[2] = 0;
		}
		else
			closest = got->second;

		if (closest[2] == 0 || (rand() % (closest[3] + closest[2])) <= closest[3])
			k = closest[0];
		else
			k = closest[1];

		closestMap[argb] = closest;
		return k;
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{		
		if (dither) 
			return dither_image(pixels, pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i, ++pixelIndex)
				qPixels[pixelIndex] = ditherFn(pPalette, nMaxColors, pixels[pixelIndex]);
		}

		return true;
	}	

	bool PnnQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);		
		
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		bool quan_sqrt = nMaxColors > BYTE_MAX;
		if (nMaxColors > 2)
			pnnquan(pixels, pPalette, nMaxColors, quan_sqrt);
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette->Entries[0] = m_transparentColor;
				pPalette->Entries[1] = Color::Black;
			}
			else {
				pPalette->Entries[0] = Color::Black;
				pPalette->Entries[1] = Color::White;
			}
		}
		
		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			dithering_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight);
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if(nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		closestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}