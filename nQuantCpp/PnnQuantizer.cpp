#pragma once
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnQuantizer.h"
#include <map>

namespace PnnQuant
{
	bool hasTransparency = false;
	ARGB m_transparentColor = Color::Transparent;
	map<ARGB, vector<short> > closestMap;
	map<ARGB, short> rightMatches;

	struct pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		int cnt = 0;
		unsigned short nn, fw, bk, tm, mtm;
	};

	inline int getIndex(const Color& c)
	{
		if (hasTransparency)
			return (c.GetA() & 0xF0) << 8 | (c.GetR() & 0xF0) << 4 | (c.GetG() & 0xF0) | (c.GetB() >> 4);
		return (c.GetR() & 0xF8) << 8 | (c.GetG() & 0xFC) << 3 | (c.GetB() >> 3);
	}

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
			double nerr, n2;

			nerr = pow((bins[i].ac - wa), 2) + pow((bins[i].rc - wr), 2) + pow((bins[i].gc - wg), 2) + pow((bins[i].bc - wb), 2);
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

	int pnnquan(const vector<ARGB>& pixels, pnnbin* bins, ColorPalette* pPalette, bool quan_sqrt)
	{
		unsigned short heap[65537] = { 0 };
		double err, n1, n2;
		int l, l2, h, b1, maxbins, extbins, res = 1;

		/* Build histogram */
		for (const auto& pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			Color c(pixel);
			int index = getIndex(c);
			auto& tb = bins[index];
			tb.ac += c.GetA();
			tb.rc += c.GetR();
			tb.gc += c.GetG();
			tb.bc += c.GetB();
			tb.cnt++;
		}

		/* Cluster nonempty bins at one end of array */
		maxbins = 0;

		for (int i = 0; i < 65536; ++i) {
			if (!bins[i].cnt)
				continue;

			double d = 1.0 / (double)bins[i].cnt;
			bins[i].ac *= d;
			bins[i].rc *= d;
			bins[i].gc *= d;
			bins[i].bc *= d;
			if (quan_sqrt)
				bins[i].cnt = sqrt(bins[i].cnt);
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
		short k = 0;
		for (int i = 0;; ++k) {
			auto alpha = rint(bins[i].ac);
			pPalette->Entries[k] = Color::MakeARGB(alpha, rint(bins[i].rc), rint(bins[i].gc), rint(bins[i].bc));
			if (hasTransparency && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);

			if (!(i = bins[i].fw))
				break;
		}

		return 0;
	}

	UINT nearestColorIndex(const ColorPalette* pPalette, const int* squares3, ARGB argb)
	{
		UINT k = 0;
		Color c(argb);

		UINT nMaxColors = pPalette->Count;
		if (hasTransparency || nMaxColors < 256) {
			auto got = rightMatches.find(argb);
			if (got == rightMatches.end()) {
				int curdist, mindist = 200000;
				const double diff = sqrt(nMaxColors);
				for (UINT i = 0; i < nMaxColors; i++) {
					Color c2(pPalette->Entries[i]);
					int adist = c2.GetA() - c.GetA();
					curdist = squares3[adist];
					if (curdist > mindist)
						continue;

					int rdist = c2.GetR() - c.GetR();
					curdist += squares3[rdist];
					if (curdist > mindist)
						continue;

					int gdist = c2.GetG() - c.GetG();
					curdist += squares3[gdist];
					if (curdist > mindist)
						continue;

					int bdist = c2.GetB() - c.GetB();
					curdist += squares3[bdist];
					if (curdist > mindist)
						continue;

					mindist = curdist;
					k = i;
				}
				rightMatches[argb] = k;
			}
			else
				k = got->second;
			return k;
		}

		vector<short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			for (; k < nMaxColors; k++) {
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

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, short* qPixels, UINT width, UINT height, bool dither)
	{
		UINT nMaxColors = pPalette->Count;
		int sqr_tbl[BYTE_MAX + BYTE_MAX + 1];

		for (int i = (-BYTE_MAX); i <= BYTE_MAX; i++)
			sqr_tbl[i + BYTE_MAX] = i * i;

		auto squares3 = &sqr_tbl[BYTE_MAX];

		UINT pixelIndex = 0;
		if (dither) {
			bool odd_scanline = false;
			short *thisrowerr, *nextrowerr;
			int j, a_pix, r_pix, g_pix, b_pix, dir, two_val;
			const byte DJ = 4;
			const byte DITHER_MAX = 20;
			const int err_len = (width + 2) * DJ;
			byte range_tbl[DJ * 256] = { 0 };
			auto range = &range_tbl[256];
			auto erowErr = make_unique<short[]>(err_len);
			auto orowErr = make_unique<short[]>(err_len);
			char dith_max_tbl[512] = { 0 };
			auto dith_max = &dith_max_tbl[256];
			auto erowerr = erowErr.get();
			auto orowerr = orowErr.get();
			short lookup[65536];

			for (int i = 0; i < 65536; i++)
				lookup[i] = -1;

			for (int i = 0; i < 256; i++) {
				range_tbl[i] = 0;
				range_tbl[i + 256] = static_cast<byte>(i);
				range_tbl[i + 512] = BYTE_MAX;
				range_tbl[i + 768] = BYTE_MAX;
			}

			for (int i = 0; i < 256; i++) {
				dith_max_tbl[i] = -DITHER_MAX;
				dith_max_tbl[i + 256] = DITHER_MAX;
			}
			for (int i = -DITHER_MAX; i <= DITHER_MAX; i++)
				dith_max_tbl[i + 256] = i;

			for (int i = 0; i < height; i++) {
				if (odd_scanline) {
					dir = -1;
					pixelIndex += (width - 1);
					thisrowerr = &orowerr[DJ];
					nextrowerr = &erowerr[width * DJ];
				}
				else {
					dir = 1;
					thisrowerr = &erowerr[DJ];
					nextrowerr = &orowerr[width * DJ];
				}
				nextrowerr[0] = nextrowerr[1] = nextrowerr[2] = nextrowerr[3] = 0;
				for (j = 0; j < width; j++) {
					Color c(pixels[pixelIndex]);

					r_pix = range[((thisrowerr[0] + 8) >> 4) + c.GetR()];
					g_pix = range[((thisrowerr[1] + 8) >> 4) + c.GetG()];
					b_pix = range[((thisrowerr[2] + 8) >> 4) + c.GetB()];
					a_pix = hasTransparency ? range[((thisrowerr[3] + 8) >> 4) + c.GetA()] : BYTE_MAX;

					ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
					Color c1(argb);
					int offset = getIndex(c1);
					if (lookup[offset] < 0)
						lookup[offset] = nearestColorIndex(pPalette, squares3, argb);
					qPixels[pixelIndex] = lookup[offset];

					Color c2(pPalette->Entries[qPixels[pixelIndex]]);
					r_pix = dith_max[r_pix - c2.GetR()];
					g_pix = dith_max[g_pix - c2.GetG()];
					b_pix = dith_max[b_pix - c2.GetB()];
					a_pix = hasTransparency ? dith_max[a_pix - c2.GetA()] : BYTE_MAX;

					two_val = r_pix * 2;
					nextrowerr[0 - DJ] = r_pix;
					r_pix += two_val;
					nextrowerr[0 + DJ] += r_pix;
					r_pix += two_val;
					nextrowerr[0] += r_pix;
					r_pix += two_val;
					thisrowerr[0 + DJ] += r_pix;

					two_val = g_pix * 2;
					nextrowerr[1 - DJ] = g_pix;
					g_pix += two_val;
					nextrowerr[1 + DJ] += g_pix;
					g_pix += two_val;
					nextrowerr[1] += g_pix;
					g_pix += two_val;
					thisrowerr[1 + DJ] += g_pix;

					two_val = b_pix * 2;
					nextrowerr[2 - DJ] = b_pix;
					b_pix += two_val;
					nextrowerr[2 + DJ] += b_pix;
					b_pix += two_val;
					nextrowerr[2] += b_pix;
					b_pix += two_val;
					thisrowerr[2 + DJ] += b_pix;

					two_val = a_pix * 2;
					nextrowerr[3 - DJ] = a_pix;
					a_pix += two_val;
					nextrowerr[3 + DJ] += a_pix;
					a_pix += two_val;
					nextrowerr[3] += a_pix;
					a_pix += two_val;
					thisrowerr[3 + DJ] += a_pix;

					thisrowerr += DJ;
					nextrowerr -= DJ;
					pixelIndex += dir;
				}
				if ((i % 2) == 1)
					pixelIndex += (width + 1);

				odd_scanline = !odd_scanline;
			}
			return true;
		}

		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++)
				qPixels[pixelIndex++] = nearestColorIndex(pPalette, squares3, pixels[pixelIndex]);
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

	bool PnnQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		if (nMaxColors > 256)
			nMaxColors = 256;
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		hasTransparency = false;
		bool r = true;
		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++) {
					Color color;
					pSource->GetPixel(x, y, &color);
					if (color.GetA() < BYTE_MAX) {
						hasTransparency = true;
						if (color.GetA() == 0)
							m_transparentColor = color.GetValue();
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
		bool quan_sqrt = nMaxColors > BYTE_MAX;
		if (nMaxColors > 2)
			pnnquan(pixels, bins.get(), pPalette, quan_sqrt);
		else {
			pPalette->Entries[0] = Color::Black;
			pPalette->Entries[1] = Color::White;
		}
		bins.reset();

		auto qPixels = make_unique<short[]>(bitmapWidth * bitmapHeight);
		quantize_image(&pixels[0], pPalette, qPixels.get(), bitmapWidth, bitmapHeight, dither);
		closestMap.clear();
		rightMatches.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}