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
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	map<ARGB, vector<short> > closestMap;

	struct pnnbin {
		double ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		int cnt = 0;
		int nn, fw, bk, tm, mtm;
	};

	inline int getARGBIndex(const Color& c)
	{
		if(hasSemiTransparency)
			return (c.GetA() & 0xF0) << 8 | (c.GetR() & 0xF0) << 4 | (c.GetG() & 0xF0) | (c.GetB() >> 4);
		if (m_transparentPixelIndex >= 0)
			return (c.GetA() & 0x80) << 8 | (c.GetR() & 0xF8) << 7 | (c.GetG() & 0xF8) << 2 | (c.GetB() >> 3);
		return (c.GetR() & 0xF8) << 8 | (c.GetG() & 0xFC) << 3 | (c.GetB() >> 3);
	}
	
	inline double sqr(double value)
    {
        return value * value;
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
			double nerr = sqr(bins[i].ac - wa) + sqr(bins[i].rc - wr) + sqr(bins[i].gc - wg) + sqr(bins[i].bc - wb);
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
		int heap[65537] = { 0 };
		double err, n1, n2;

		/* Build histogram */
		for (const auto& pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?
			Color c(pixel);
			int index = getARGBIndex(c);
			auto& tb = bins[index];
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
			if (m_transparentPixelIndex >= 0 && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);

			if (!(i = bins[i].fw))
				break;
		}

		return 0;
	}

	UINT nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		UINT k = 0;
		Color c(argb);

		UINT curdist, mindist = SHORT_MAX;
		for (short i = 0; i < nMaxColors; i++) {
			Color c2(pPalette->Entries[i]);
			int adist = abs(c2.GetA() - c.GetA());
			curdist = adist;
			if (curdist > mindist)
				continue;

			int rdist = abs(c2.GetR() - c.GetR());
			curdist += rdist;
			if (curdist > mindist)
				continue;

			int gdist = abs(c2.GetG() - c.GetG());
			curdist += gdist;
			if (curdist > mindist)
				continue;

			int bdist = abs(c2.GetB() - c.GetB());
			curdist += bdist;
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		return k;
	}

	UINT closestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		UINT k = 0;
		Color c(argb);
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

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		UINT pixelIndex = 0;
		if (dither) {
			bool odd_scanline = false;
			short *row0, *row1;
			int a_pix, r_pix, g_pix, b_pix, dir, k;
			const int DJ = 4;
			const int DITHER_MAX = 20;
			const int err_len = (width + 2) * DJ;
			byte clamp[DJ * 256] = { 0 };
			auto erowErr = make_unique<short[]>(err_len);
			auto orowErr = make_unique<short[]>(err_len);
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

			for (UINT i = 0; i < height; i++) {
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
				for (UINT j = 0; j < width; j++) {
					Color c(pixels[pixelIndex]);

					r_pix = clamp[((row0[0] + 0x1008) >> 4) + c.GetR()];
					g_pix = clamp[((row0[1] + 0x1008) >> 4) + c.GetG()];
					b_pix = clamp[((row0[2] + 0x1008) >> 4) + c.GetB()];
					a_pix = clamp[((row0[3] + 0x1008) >> 4) + c.GetA()];

					ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
					Color c1(argb);
					int offset = getARGBIndex(c1);
					if (!lookup[offset])
						lookup[offset] = nearestColorIndex(pPalette, nMaxColors, c.GetA() ? argb : pixels[pixelIndex]) + 1;
					qPixels[pixelIndex] = lookup[offset] - 1;

					Color c2(pPalette->Entries[qPixels[pixelIndex]]);

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

		UINT(*fcnPtr)(const ColorPalette*, const UINT nMaxColors, const ARGB) = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		for (int j = 0; j < height; j++) {
			for (int i = 0; i < width; i++)
				qPixels[pixelIndex++] = (*fcnPtr)(pPalette, nMaxColors, pixels[pixelIndex]);
		}

		return true;
	}

	bool quantize_image(const ARGB* pixels, short* qPixels, UINT width, UINT height)
	{
		UINT nMaxColors = 65536;

		UINT pixelIndex = 0;
		bool odd_scanline = false;
		short *row0, *row1;
		int a_pix, r_pix, g_pix, b_pix, dir, k;
		const int DJ = 4;
		const int DITHER_MAX = 20;
		const int err_len = (width + 2) * DJ;
		byte clamp[DJ * 256] = { 0 };
		auto erowErr = make_unique<short[]>(err_len);
		auto orowErr = make_unique<short[]>(err_len);
		char limtb[512] = { 0 };
		auto lim = &limtb[256];
		auto erowerr = erowErr.get();
		auto orowerr = orowErr.get();
		int lookup[65536] = { 0 };

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
			for (UINT j = 0; j < width; j++) {
				Color c(pixels[pixelIndex]);

				r_pix = clamp[((row0[0] + 0x1008) >> 4) + c.GetR()];
				g_pix = clamp[((row0[1] + 0x1008) >> 4) + c.GetG()];
				b_pix = clamp[((row0[2] + 0x1008) >> 4) + c.GetB()];
				a_pix = clamp[((row0[3] + 0x1008) >> 4) + c.GetA()];

				ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
				Color c1(argb);
				int offset = getARGBIndex(c1);
				if (!lookup[offset]) {
					auto argb1 = Color::MakeARGB(BYTE_MAX, (c1.GetR() & 0xF8), (c1.GetG() & 0xFC), (c1.GetB() & 0xF8));
					if (hasSemiTransparency)
						argb1 = Color::MakeARGB((c1.GetA() & 0xF0), (c1.GetR() & 0xF0), (c1.GetG() & 0xF0), (c1.GetB() & 0xF0));
					else if (m_transparentPixelIndex >= 0)
						argb1 = Color::MakeARGB((c1.GetA() < BYTE_MAX) ? 0 : BYTE_MAX, (c1.GetR() & 0xF8), (c1.GetG() & 0xF8), (c1.GetB() & 0xF8));
					lookup[offset] = argb1 + 1;
				}
				qPixels[pixelIndex] = static_cast<short>(offset);

				Color c2(static_cast<ARGB>(lookup[offset] - 1));

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

	bool ProcessImagePixels(Bitmap* pDest, const short* qPixels)
	{
		UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
		if (bpp < 16)
			return false;

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

		// Second loop: fill indexed bitmap
		for (UINT y = 0; y < h; y++) {	// For each row...
			for (UINT x = 0; x < w * 2;) {
				auto argb = static_cast<unsigned short>(qPixels[pixelIndex++]);
				pRowDest[x++] = static_cast<byte>(argb & 0xFF);
				pRowDest[x++] = static_cast<byte>(argb >> 8);				
			}
			pRowDest += strideDest;
		}

		status = pDest->UnlockBits(&targetData);
		return pDest->GetLastStatus() == Ok;
	}

	bool PnnQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		hasSemiTransparency = false;
		m_transparentPixelIndex = -1;
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
							m_transparentPixelIndex = pixelIndex;
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
						hasSemiTransparency = true;
						if (pixelAlpha == 0) {
							m_transparentPixelIndex = pixelIndex;
							m_transparentColor = argb;
						}
					}
					pixels[pixelIndex++] = argb;
				}

				pRowSource += strideSource;
			}

			pSource->UnlockBits(&data);
		}

		if (nMaxColors > 256) {
			hasSemiTransparency = false;
			auto qPixels = make_unique<short[]>(bitmapWidth * bitmapHeight);
			quantize_image(pixels.data(), qPixels.get(), bitmapWidth, bitmapHeight);
			return ProcessImagePixels(pDest, qPixels.get());
		}
		
		auto pPaletteBytes = make_unique<byte[]>(pDest->GetPaletteSize());
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

		auto qPixels = make_unique<short[]>(bitmapWidth * bitmapHeight);
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