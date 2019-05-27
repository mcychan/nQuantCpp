#pragma once
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2019 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnLABQuantizer.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include <unordered_map>

namespace PnnLABQuant
{
	double PR = .2126, PG = .7152, PB = .0722;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
	unordered_map<ARGB, vector<double> > closestMap;

	struct pnnbin {
		double ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
		int cnt = 0;
		int nn = 0, fw = 0, bk = 0, tm = 0, mtm = 0;
	};

	void getLab(const Color& c, CIELABConvertor::Lab& lab1)
	{
		auto got = pixelMap.find(c.GetValue());
		if (got == pixelMap.end()) {
			CIELABConvertor::RGB2LAB(c, lab1);
			pixelMap[c.GetValue()] = lab1;
		}
		else
			lab1 = got->second;
	}

	void find_nn(pnnbin* bins, int idx)
	{
		int i, nn = 0;
		double err = 1e100;

		auto& bin1 = bins[idx];
		auto n1 = bin1.cnt;
		auto wa = bin1.ac;
		auto wL = bin1.Lc;
		auto wA = bin1.Ac;
		auto wB = bin1.Bc;
		for (i = bin1.fw; i; i = bins[i].fw) {
			double nerr, n2;

			nerr = sqr(bins[i].ac - wa) + sqr(bins[i].Lc - wL) + sqr(bins[i].Ac - wA) + sqr(bins[i].Bc - wB);
			n2 = bins[i].cnt;
			nerr *= (n1 * n2) / (n1 + n2);
			if (nerr >= err)
				continue;
			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	int PnnLABQuantizer::pnnquan(const vector<ARGB>& pixels, ColorPalette* pPalette, UINT nMaxColors)
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

			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
			auto& tb = bins[index];
			tb.ac += c.GetA();
			tb.Lc += lab1.L;
			tb.Ac += lab1.A;
			tb.Bc += lab1.B;
			tb.cnt++;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < 65536; ++i) {
			if (!bins[i].cnt)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].Lc *= d;
			bins[i].Ac *= d;
			bins[i].Bc *= d;
			bins[maxbins++] = bins[i];
		}

		for (int i = 0; i < maxbins - 1; i++) {
			bins[i].fw = i + 1;
			bins[i + 1].bk = i;
		}
		// !!! Already zeroed out by calloc()
		//	bins[0].bk = bins[i].fw = 0;

		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; i++) {
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
			tb.ac = d * (n1 * tb.ac + n2 * nb.ac);
			tb.Lc = d * (n1 * tb.Lc + n2 * nb.Lc);
			tb.Ac = d * (n1 * tb.Ac + n2 * nb.Ac);
			tb.Bc = d * (n1 * tb.Bc + n2 * nb.Bc);
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = 0xFFFF;
		}

		/* Fill palette */
		short k = 0;

		for (int i = 0;; k++) {			
			CIELABConvertor::Lab lab1;
			lab1.alpha = rint(bins[i].ac);
			lab1.L = bins[i].Lc, lab1.A = bins[i].Ac, lab1.B = bins[i].Bc;
			pPalette->Entries[k] = CIELABConvertor::LAB2RGB(lab1);
			if (m_transparentPixelIndex >= 0 && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);

			if (!(i = bins[i].fw))
				break;
		}

		return 0;
	}

	short nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		short k = 0;
		Color c(argb);

		double mindist = INT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (short i = 0; i < nMaxColors; i++) {
			Color c2(pPalette->Entries[i]);
			if (nMaxColors > 32) {
				double curdist = sqr(c2.GetA() - c.GetA());
				if (curdist > mindist)
					continue;

				curdist += PR * sqr(c2.GetR() - c.GetR());
				if (curdist > mindist)
					continue;

				curdist += PG * sqr(c2.GetG() - c.GetG());
				if (curdist > mindist)
					continue;

				curdist += PB * sqr(c2.GetB() - c.GetB());
				if (curdist > mindist)
					continue;

				mindist = curdist;
			}
			else {
				getLab(c2, lab2);

				double curdist = sqr(c2.GetA() - c.GetA());
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.B - lab1.B);
				if (curdist > mindist)
					continue;

				mindist = curdist;
			}
			
			k = i;
		}
		return k;
	}
	
	short closestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		short k = 0;
		Color c(argb);
		vector<double> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			CIELABConvertor::Lab lab1, lab2;
			getLab(c, lab1);
			for (; k < nMaxColors; k++) {
				Color c2(pPalette->Entries[k]);
				getLab(c2, lab2);
				closest[4] = sqr(lab2.alpha - lab1.alpha) + CIELABConvertor::CIEDE2000(lab2, lab1);
				//closest[4] = abs(lab2.alpha - lab1.alpha) + abs(lab2.L - lab1.L) + abs(lab2.A - lab1.A) + abs(lab2.B - lab1.B);
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

		if (closest[2] == 0 || (rand() % (int) ceil(closest[3] + closest[2])) <= closest[3])
			k = closest[0];
		else
			k = closest[1];

		closestMap[argb] = closest;
		return k;
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither) 
			return dither_image(pixels, pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i)
				qPixels[pixelIndex++] = ditherFn(pPalette, nMaxColors, pixels[pixelIndex]);
		}
		return true;
	}

	bool PnnLABQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);
		
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors > 2)
			pnnquan(pixels, pPalette, nMaxColors);
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
		
		auto qPixels = make_unique<short[]>(pixels.size());
		if (nMaxColors > 256) {
			hasSemiTransparency = false;
			dithering_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight);
			return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex);
		}
		if (hasSemiTransparency)
			PR = PG = PB = 1;
		quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);
		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		pixelMap.clear();
		closestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}
