#pragma once
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnLABQuantizer.h"
#include "CIELABConvertor.h"
#include <map>

namespace PnnLABQuant
{
	bool hasTransparency = false, hasSemiTransparency = false;
	ARGB m_transparentColor = Color::Transparent;
	map<ARGB, CIELABConvertor::Lab> pixelMap;	
	map<ARGB, vector<double> > closestMap;

	inline int getARGBIndex(const Color& c)
	{
		if (hasSemiTransparency)
			return (c.GetA() & 0xF0) << 8 | (c.GetR() & 0xF0) << 4 | (c.GetG() & 0xF0) | (c.GetB() >> 4);
		if (hasTransparency)
			return (c.GetA() & 0x80) << 8 | (c.GetR() & 0xF8) << 7 | (c.GetG() & 0xF8) << 2 | (c.GetB() >> 3);
		return (c.GetR() & 0xF8) << 8 | (c.GetG() & 0xFC) << 3 | (c.GetB() >> 3);
	}

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

			nerr = pow((bins[i].ac - wa), 2) + pow((bins[i].Lc - wL), 2) + pow((bins[i].Ac - wA), 2) + pow((bins[i].Bc - wB), 2);
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

	int PnnLABQuantizer::pnnquan(const vector<ARGB>& pixels, pnnbin* bins, ColorPalette* pPalette)
	{
		unsigned short heap[65537] = { 0 };
		double err, n1, n2;
		int l, l2, h, b1, maxbins, extbins, res = 1;

		/* Build histogram */
		for (const auto& pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?			
			Color c(pixel);
			int index = getARGBIndex(c);

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
		maxbins = 0;

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

		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; i++) {
			find_nn(bins, i);
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
		extbins = maxbins - pPalette->Count;
		for (int i = 0; i < extbins; ) {
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
					find_nn(bins, b1);
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
		if (hasTransparency)
			++k;
		for (int i = 0;; k++) {			
			CIELABConvertor::Lab lab1;
			lab1.alpha = rint(bins[i].ac);
			lab1.L = bins[i].Lc, lab1.A = bins[i].Ac, lab1.B = bins[i].Bc;
			pPalette->Entries[k] = CIELABConvertor::LAB2RGB(lab1);
			if (hasTransparency && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);

			if (!(i = bins[i].fw))
				break;
		}

		return 0;
	}

	short nearestColorIndex(const ColorPalette* pPalette, ARGB argb)
	{
		short k = 0;
		Color c(argb);

		unsigned short index = getARGBIndex(c);
		double mindist = SHORT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = 0; i < pPalette->Count; i++) {
			Color c2(pPalette->Entries[i]);
			getLab(c2, lab2);

			double curdist = pow(c2.GetA() - c.GetA(), 2.0);
			if (curdist > mindist)
				continue;

			if (pPalette->Count < 256) {
				double deltaL_prime_div_k_L_S_L = CIELABConvertor::L_prime_div_k_L_S_L(lab1, lab2);
				curdist += pow(deltaL_prime_div_k_L_S_L, 2.0);
				if (curdist > mindist)
					continue;

				double a1Prime, a2Prime, CPrime1, CPrime2;
				double deltaC_prime_div_k_L_S_L = CIELABConvertor::C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += pow(deltaC_prime_div_k_L_S_L, 2.0);
				if (curdist > mindist)
					continue;

				double barCPrime, barhPrime;
				double deltaH_prime_div_k_L_S_L = CIELABConvertor::H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist += pow(deltaH_prime_div_k_L_S_L, 2.0);
				if (curdist > mindist)
					continue;

				curdist += CIELABConvertor::R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
			}
			else {
				curdist += pow(lab2.L - lab1.L, 2.0);
				if (curdist > mindist)
					continue;

				curdist += pow(lab2.A - lab1.A, 2.0);
				if (curdist > mindist)
					continue;

				curdist += pow(lab2.B - lab1.B, 2.0);
				if (curdist > mindist)
					continue;
			}

			mindist = curdist;
			k = i;
		}
		return k;
	}
	
	short closestColorIndex(const ColorPalette* pPalette, const ARGB argb)
	{
		short k = 0;
		Color c(argb);
		vector<double> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			CIELABConvertor::Lab lab1, lab2;
			getLab(c, lab1);
			UINT nMaxColors = pPalette->Count;
			for (; k < nMaxColors; k++) {
				Color c2(pPalette->Entries[k]);
				getLab(c2, lab2);
				closest[4] = pow(lab2.alpha - lab1.alpha, 2) + CIELABConvertor::CIEDE2000(lab2, lab1);
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

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, short* qPixels, UINT nMaxColors, int width, int height, bool dither)
	{
		UINT pixelIndex = 0;

		if (dither) {
			bool odd_scanline = false;
			short *row0, *row1;
			int j, a_pix, r_pix, g_pix, b_pix, dir, k;
			const int DJ = 4;
			const int DITHER_MAX = 20;
			const int err_len = (width + 2) * DJ;
			byte clamp[DJ * 256] = { 0 };
			unique_ptr<short[]> erowErr = make_unique<short[]>(err_len);
			unique_ptr<short[]> orowErr = make_unique<short[]>(err_len);
			char limtb[512] = { 0 };
			auto lim = &limtb[256];
			auto erowerr = erowErr.get();
			auto orowerr = orowErr.get();
			short lookup[65536] = { 0 };

			for (int i = 0; i < 256; i++) {
				clamp[i] = 0;
				clamp[i + 256] = static_cast<byte>(i);
				clamp[i + 512] = BYTE_MAX;
				clamp[i + 768] = BYTE_MAX;

				limtb[i] = -DITHER_MAX;
				limtb[i + 256] = DITHER_MAX;
			}
			for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
				limtb[i + 256] = i;

			for (int i = 0; i < height; i++) {
				if (odd_scanline) {
					dir = -1;
					pixelIndex += (width - 1);
					row0 = &orowerr[DJ];
					row1 = &erowerr[width * DJ];
				}
				else {
					dir = 1;
					row0 = &erowerr[DJ];
					row1 = &orowerr[width * DJ];
				}
				row1[0] = row1[1] = row1[2] = row1[3] = 0;
				for (j = 0; j < width; j++) {
					Color c(pixels[pixelIndex]);

					r_pix = clamp[((row0[0] + 0x1008) >> 4) + c.GetR()];
					g_pix = clamp[((row0[1] + 0x1008) >> 4) + c.GetG()];
					b_pix = clamp[((row0[2] + 0x1008) >> 4) + c.GetB()];
					a_pix = clamp[((row0[3] + 0x1008) >> 4) + c.GetA()];

					ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
					Color c1(argb);
					int offset = getARGBIndex(c1);
					if (!lookup[offset])
						lookup[offset] = nearestColorIndex(pPalette, argb) + 1;
					qPixels[pixelIndex] = lookup[offset] - 1;

					Color c2(pPalette->Entries[qPixels[pixelIndex]]);
					if (c2.GetA() < BYTE_MAX && c.GetA() == BYTE_MAX) {
						lookup[offset] = nearestColorIndex(pPalette, pixels[pixelIndex]) + 1;
						qPixels[pixelIndex] = lookup[offset] - 1;
					}

					r_pix = lim[r_pix - c2.GetR()];
					g_pix = lim[g_pix - c2.GetG()];
					b_pix = lim[b_pix - c2.GetB()];
					a_pix = lim[a_pix - c2.GetA()];

					k = r_pix * 2;
					row1[0 - DJ] = r_pix;
					row1[0 + DJ] += (r_pix += k);
					row1[0] += (r_pix += k);
					row0[0 + DJ] += (r_pix += k);

					k = g_pix * 2;
					row1[1 - DJ] = g_pix;
					row1[1 + DJ] += (g_pix += k);
					row1[1] += (g_pix += k);
					row0[1 + DJ] += (g_pix += k);

					k = b_pix * 2;
					row1[2 - DJ] = b_pix;
					row1[2 + DJ] += (b_pix += k);
					row1[2] += (b_pix += k);
					row0[2 + DJ] += (b_pix += k);

					k = a_pix * 2;
					row1[3 - DJ] = a_pix;
					row1[3 + DJ] += (a_pix += k);
					row1[3] += (a_pix += k);
					row0[3 + DJ] += (a_pix += k);

					row0 += DJ;
					row1 -= DJ;
					pixelIndex += dir;
				}
				if ((i % 2) == 1)
					pixelIndex += (width + 1);

				odd_scanline = !odd_scanline;
			}
			return true;
		}

		short(*fcnPtr)(const ColorPalette*, const ARGB) = (hasTransparency || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++)
				qPixels[pixelIndex++] = (*fcnPtr)(pPalette, pixels[pixelIndex]);
		}
		return true;
	}

	bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, const short* qPixels)
	{
		pDest->SetPalette(pPalette);

		BitmapData targetData;
		UINT w = pDest->GetWidth();
		UINT h = pDest->GetHeight();

		Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
		if (status != Ok) {
			AfxMessageBox(_T("Cannot write image"));
			return false;
		}

		int pixelIndex = 0;

		auto pRowDest = (byte*)targetData.Scan0;
		UINT strideDest;

		// Compensate for possible negative stride
		if (targetData.Stride > 0)
			strideDest = targetData.Stride;
		else {
			pRowDest += h * targetData.Stride;
			strideDest = -targetData.Stride;
		}

		UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
		// Second loop: fill indexed bitmap
		for (UINT y = 0; y < h; y++) {	// For each row...
			for (UINT x = 0; x < w; x++) {	// ...for each pixel...
				byte nibbles = 0;
				byte index = static_cast<byte>(qPixels[pixelIndex++]);

				switch (bpp)
				{
				case 8:
					pRowDest[x] = index;
					break;
				case 4:
					// First pixel is the high nibble. From and To indices are 0..16
					nibbles = pRowDest[x / 2];
					if ((x & 1) == 0) {
						nibbles &= 0x0F;
						nibbles |= (byte)(index << 4);
					}
					else {
						nibbles &= 0xF0;
						nibbles |= index;
					}

					pRowDest[x / 2] = nibbles;
					break;
				case 1:
					// First pixel is MSB. From and To are 0 or 1.
					int pos = x / 8;
					byte mask = (byte)(128 >> (x & 7));
					if (index == 0)
						pRowDest[pos] &= (byte)~mask;
					else
						pRowDest[pos] |= mask;
					break;
				}
			}

			pRowDest += strideDest;
		}

		status = pDest->UnlockBits(&targetData);
		return pDest->GetLastStatus() == Ok;
	}

	bool PnnLABQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		if (nMaxColors > 256)
			nMaxColors = 256;
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		hasTransparency = hasSemiTransparency = false;
		bool r = true;
		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++) {
					Color color;
					pSource->GetPixel(x, y, &color);
					if (color.GetA() < BYTE_MAX) {
						hasSemiTransparency = true;
						if (color.GetA() == 0) {
							hasTransparency = true;
							m_transparentColor = color.GetValue();
						}
					}
					pixels[pixelIndex++] = color.GetValue();
				}
			}
		}

		// Lock bits on 3x8 source bitmap
		else {
			BitmapData data;
			Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, pSource->GetPixelFormat(), &data);
			if (status != Ok)
				return false;

			auto pRowSource = (byte*)data.Scan0;
			UINT strideSource;

			if (data.Stride > 0) strideSource = data.Stride;
			else
			{
				// Compensate for possible negative stride
				// (not needed for first loop, but we have to do it
				// for second loop anyway)
				pRowSource += bitmapHeight * data.Stride;
				strideSource = -data.Stride;
			}

			int pixelIndex = 0;

			// First loop: gather color information
			for (UINT y = 0; y < bitmapHeight; y++) {	// For each row...
				auto pPixelSource = pRowSource;

				for (UINT x = 0; x < bitmapWidth; x++) {	// ...for each pixel...
					byte pixelBlue = *pPixelSource++;
					byte pixelGreen = *pPixelSource++;
					byte pixelRed = *pPixelSource++;
					byte pixelAlpha = bitDepth < 32 ? BYTE_MAX : *pPixelSource++;

					auto argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
					if (pixelAlpha < BYTE_MAX) {
						hasTransparency = true;
						if (pixelAlpha == 0)
							m_transparentColor = argb;
					}
					pixels[pixelIndex++] = argb;
				}

				pRowSource += strideSource;
			}

			pSource->UnlockBits(&data);
		}

		auto bins = make_unique<pnnbin[]>(65536);
		auto pPaletteBytes = make_unique<byte[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		if (nMaxColors > 2)
			pnnquan(pixels, bins.get(), pPalette);
		else {
			if (hasTransparency) {
				pPalette->Entries[0] = Color::Transparent;
				pPalette->Entries[1] = Color::Black;
			}
			else {
				pPalette->Entries[0] = Color::Black;
				pPalette->Entries[1] = Color::White;
			}
		}
		bins.reset();

		auto qPixels = make_unique<short[]>(bitmapWidth * bitmapHeight);
		quantize_image(pixels.data(), pPalette, qPixels.get(), nMaxColors, bitmapWidth, bitmapHeight, dither);
		pixelMap.clear();
		closestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}
