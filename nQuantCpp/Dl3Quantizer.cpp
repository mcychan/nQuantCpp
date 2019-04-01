#pragma once
/*
 * DL3 Quantization
 * ================
 *
 * File: dl3quant.c
 * Author: Dennis Lee   E-mail: denlee@ecf.utoronto.ca
 *
 * Copyright (C) 1993-1997 Dennis Lee
 *
 * C implementation of DL3 Quantization.
 * DL3 Quantization is a 2-pass color quantizer that uses an
 * exhaustive search technique to minimize error introduced at
 * each step during palette reduction.
 *
 * I believe DL3 Quant offers the highest quality of all existing
 * color quantizers.  It is truly 'optimal' except for a few provisos.
 * These provisos and other information about DL3 Quant can be found
 * in DLQUANT.TXT, which is included in this distribution.
 *
 *
 * NOTES
 * =====
 *
 * The dithering code is based on code from the IJG's jpeg library.
 *
 * DL3 Quantization can take a long time to reduce a palette.
 * Times can range from seconds to minutes or even hours depending on
 * the image and the computer used.  This eliminates DL3 Quant for
 * typical usage unless the user has a very fast computer and/or has a
 * lot of patience.  However, the reward is a quantized image that is
 * the best currently possible.  The number of colors in the source image,
 * not the image size, determines the time required to quantize it.
 *
 * This source code may be freely copied, modified, and redistributed,
 * provided this copyright notice is attached.
 * Compiled versions of this code, modified or not, are free for
 * personal use.  Compiled versions used in distributed software
 * is also free, but a notification must be sent to the author.
 * An e-mail to denlee@ecf.utoronto.ca will do.
 *
 */

#include "stdafx.h"
#include "Dl3Quantizer.h"
#include "bitmapUtilities.h"
#include <map>

namespace Dl3Quant
{
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	map<ARGB, vector<short> > closestMap;
	map<ARGB, UINT> rightMatches;

	struct CUBE3 {
		UINT a, r, g, b;
		UINT pixel_count = 0;
		UINT err;
		int cc;
		byte aa, rr, gg, bb;
	};

	void setrgb(CUBE3& rec)
	{
		int v = rec.pixel_count, v2 = v >> 1;
		rec.aa = (rec.a + v2) / v;
		rec.rr = (rec.r + v2) / v;
		rec.gg = (rec.g + v2) / v;
		rec.bb = (rec.b + v2) / v;
	}
	
	UINT calc_err(CUBE3* rgb_table3, int* squares3, int c1, int c2)
	{
		UINT dist1, dist2, P1, P2, P3;
		int A1, R1, G1, B1, A2, R2, G2, B2, A3, R3, G3, B3;

		P1 = rgb_table3[c1].pixel_count;
		P2 = rgb_table3[c2].pixel_count;
		P3 = P1 + P2;

		A3 = (rgb_table3[c1].a + rgb_table3[c2].a + (P3 >> 1)) / P3;
		R3 = (rgb_table3[c1].r + rgb_table3[c2].r + (P3 >> 1)) / P3;
		G3 = (rgb_table3[c1].g + rgb_table3[c2].g + (P3 >> 1)) / P3;
		B3 = (rgb_table3[c1].b + rgb_table3[c2].b + (P3 >> 1)) / P3;

		A1 = rgb_table3[c1].aa;
		R1 = rgb_table3[c1].rr;
		G1 = rgb_table3[c1].gg;
		B1 = rgb_table3[c1].bb;

		A2 = rgb_table3[c2].aa;
		R2 = rgb_table3[c2].rr;
		G2 = rgb_table3[c2].gg;
		B2 = rgb_table3[c2].bb;

		dist1 = squares3[A3 - A1] + squares3[R3 - R1] + squares3[G3 - G1] + squares3[B3 - B1];
		dist1 = (UINT)(sqrt(dist1) * P1);

		dist2 = squares3[A2 - A3] + squares3[R2 - R3] + squares3[G2 - G3] + squares3[B2 - B3];
		dist2 = (UINT)(sqrt(dist2) * P2);

		return (dist1 + dist2);
	}

	void build_table3(CUBE3* rgb_table3, ARGB argb)
	{
		Color c(argb);
		int index = getARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex);

		rgb_table3[index].a += c.GetA();
		rgb_table3[index].r += c.GetR();
		rgb_table3[index].g += c.GetG();
		rgb_table3[index].b += c.GetB();
		rgb_table3[index].pixel_count++;
	}

	UINT build_table3(CUBE3* rgb_table3, vector<ARGB>& pixels)
	{
		for (const auto & pixel : pixels)
			build_table3(rgb_table3, pixel);

		UINT tot_colors = 0;
		for (int i = 0; i < 65536; i++) {			
			if (rgb_table3[i].pixel_count) {
				setrgb(rgb_table3[i]);
				rgb_table3[tot_colors++] = rgb_table3[i];
			}
		}
		return tot_colors;
	}

	void recount_next(CUBE3* rgb_table3, int* squares3, UINT tot_colors, int i)
	{
		int c2 = 0;
		UINT err, cur_err;

		err = ~0L;
		for (int j = i + 1; j < tot_colors; j++) {
			cur_err = calc_err(rgb_table3, squares3, i, j);
			if (cur_err < err) {
				err = cur_err;
				c2 = j;
			}
		}
		rgb_table3[i].err = err;
		rgb_table3[i].cc = c2;
	}

	void recount_dist(CUBE3* rgb_table3, int* squares3, UINT tot_colors, int c1)
	{
		int i;
		UINT cur_err;

		recount_next(rgb_table3, squares3, tot_colors, c1);
		for (i = 0; i < c1; i++) {
            if (rgb_table3[i].cc == c1)
				recount_next(rgb_table3, squares3, tot_colors, i);
			else {
				cur_err = calc_err(rgb_table3, squares3, i, c1);
				if (cur_err < rgb_table3[i].err) {
					rgb_table3[i].err = cur_err;
					rgb_table3[i].cc = c1;
				}
			}
		}
	}

	int reduce_table3(CUBE3* rgb_table3, int* squares3, UINT tot_colors, UINT num_colors)
	{
		int c1=0, c2=0, grand_total, bailout = FALSE;
		UINT err;

		int i = 0;
		for (; i < (tot_colors - 1); i++)
			recount_next(rgb_table3, squares3, tot_colors, i);

		rgb_table3[i].err = ~0L;
		rgb_table3[i].cc = tot_colors;

		grand_total = tot_colors-num_colors;
		while (tot_colors > num_colors) {
			err = ~0L;
			for (int i = 0; i < tot_colors; i++) {
                  if (rgb_table3[i].err < err)
                  {
                        err = rgb_table3[i].err;
                        c1 = i;
                  }
            }
            c2 = rgb_table3[c1].cc;
			rgb_table3[c2].a += rgb_table3[c1].a;
            rgb_table3[c2].r += rgb_table3[c1].r;
            rgb_table3[c2].g += rgb_table3[c1].g;
            rgb_table3[c2].b += rgb_table3[c1].b;
            rgb_table3[c2].pixel_count += rgb_table3[c1].pixel_count;
            setrgb(rgb_table3[c2]);
            tot_colors--;

            rgb_table3[c1] = rgb_table3[tot_colors];
            rgb_table3[tot_colors-1].err = ~0L;
            rgb_table3[tot_colors-1].cc = tot_colors;

            for (i = 0; i < c1; i++) {
				if (rgb_table3[i].cc == tot_colors)
					rgb_table3[i].cc = c1;
            }

            for (i = c1 + 1; i < tot_colors; i++) {
				if (rgb_table3[i].cc == tot_colors)
					recount_next(rgb_table3, squares3, tot_colors, i);
            }

            recount_dist(rgb_table3, squares3, tot_colors, c1);
            if (c2 != tot_colors)
				recount_dist(rgb_table3, squares3, tot_colors, c2);
		}

		return bailout;
	}
	
	short nearestColorIndex(const CUBE3* rgb_table3, const int* squares3, const UINT nMaxColors, const ARGB argb)
	{
		short k = 0;
		Color c(argb);

		auto got = rightMatches.find(argb);
		if (got == rightMatches.end()) {
			int curdist, mindist = 200000;
			for (UINT i = 0; i < nMaxColors; i++) {
				curdist = squares3[rgb_table3[i].aa - c.GetA()];
				if (curdist > mindist)
					continue;
				curdist += squares3[rgb_table3[i].rr - c.GetR()];
				if (curdist > mindist)
					continue;
				curdist += squares3[rgb_table3[i].gg - c.GetG()];
				if (curdist > mindist)
					continue;
				curdist += squares3[rgb_table3[i].bb - c.GetB()];
				if (curdist > mindist)
					continue;
				if (curdist < mindist) {
					mindist = curdist;					
					k = i;
				}
			}
			rightMatches[argb] = k;
		}
		else
			k = got->second;

		return k;
	}
		
	short closestColorIndex(const CUBE3* rgb_table3, const int* squares3, const UINT nMaxColors, const ARGB argb)
	{
		short k = 0;
		Color c(argb);
		vector<short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			for (; k < nMaxColors; k++) {
				closest[4] = abs(c.GetA() - rgb_table3[k].aa) + abs(c.GetR() - rgb_table3[k].rr) + abs(c.GetG() - rgb_table3[k].gg) + abs(c.GetB() - rgb_table3[k].bb);
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

			if (closest[3] == 100000000)
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
	
	bool quantize_image3(const ARGB* pixels, CUBE3* rgb_table3, const UINT nMaxColors, short* qPixels, int* squares3, const UINT width, const UINT height, const bool dither)
	{			
		if (dither) {
			bool odd_scanline = false;
			short *thisrowerr, *nextrowerr;
			int j, a_pix, r_pix, g_pix, b_pix, dir, two_val;
			const byte DJ = 4;
			const byte DITHER_MAX = 20;
			const int err_len = (width + 2) * DJ;
			byte range_tbl[DJ * 256] = { 0 };
			auto range = &range_tbl[256];
			unique_ptr<short[]> erowErr = make_unique<short[]>(err_len);
			unique_ptr<short[]> orowErr = make_unique<short[]>(err_len);
			char dith_max_tbl[512] = { 0 };
			auto dith_max = &dith_max_tbl[256];
			auto erowerr = erowErr.get();
			auto orowerr = orowErr.get();

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

			UINT pixelIndex = 0;
			for (int i = 0; i < height; i++) {
				if (odd_scanline) {
					dir = -1;
					pixelIndex += (width - 1);
					thisrowerr = orowerr + DJ;
					nextrowerr = erowerr + width * DJ;
				}
				else {
					dir = 1;
					thisrowerr = erowerr + DJ;
					nextrowerr = orowerr + width * DJ;
				}
				nextrowerr[0] = nextrowerr[1] = nextrowerr[2] = nextrowerr[3] = 0;
				for (j = 0; j < width; j++) {
					Color c(pixels[pixelIndex]);

					a_pix = range[((thisrowerr[0] + 8) >> 4) + c.GetA()];
					r_pix = range[((thisrowerr[1] + 8) >> 4) + c.GetR()];
					g_pix = range[((thisrowerr[2] + 8) >> 4) + c.GetG()];
					b_pix = range[((thisrowerr[3] + 8) >> 4) + c.GetB()];

					ARGB argb = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
					qPixels[pixelIndex] = nearestColorIndex(rgb_table3, squares3, nMaxColors, argb);

					Color c2(qPixels[pixelIndex]);
					a_pix = dith_max[a_pix - c2.GetA()];
					r_pix = dith_max[r_pix - c2.GetR()];
					g_pix = dith_max[g_pix - c2.GetG()];
					b_pix = dith_max[b_pix - c2.GetB()];

					two_val = a_pix * 2;
					nextrowerr[0 - DJ] = a_pix;
					a_pix += two_val;
					nextrowerr[0 + DJ] += a_pix;
					a_pix += two_val;
					nextrowerr[0] += a_pix;
					a_pix += two_val;
					thisrowerr[0 + DJ] += a_pix;

					two_val = r_pix * 2;
					nextrowerr[1 - DJ] = r_pix;
					r_pix += two_val;
					nextrowerr[1 + DJ] += r_pix;
					r_pix += two_val;
					nextrowerr[1] += r_pix;
					r_pix += two_val;
					thisrowerr[1 + DJ] += r_pix;

					two_val = g_pix * 2;
					nextrowerr[2 - DJ] = g_pix;
					g_pix += two_val;
					nextrowerr[2 + DJ] += g_pix;
					g_pix += two_val;
					nextrowerr[2] += g_pix;
					g_pix += two_val;
					thisrowerr[2 + DJ] += g_pix;

					two_val = b_pix * 2;
					nextrowerr[3 - DJ] = b_pix;
					b_pix += two_val;
					nextrowerr[3 + DJ] += b_pix;
					b_pix += two_val;
					nextrowerr[3] += b_pix;
					b_pix += two_val;
					thisrowerr[3 + DJ] += b_pix;

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

		short(*fcnPtr)(const CUBE3* rgb_table3, const int* squares3, const UINT nMaxColors, const ARGB) = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		for (int i = 0; i < (width * height); ++i)
			qPixels[i] = (*fcnPtr)(rgb_table3, squares3, nMaxColors, pixels[i]);
		return true;
	}

	void GetQuantizedPalette(ColorPalette* pPalette, const CUBE3* rgb_table3)
	{
		short k = 0;
		for (UINT paletteIndex = 0; paletteIndex < pPalette->Count; paletteIndex++) {
			UINT sum = rgb_table3[paletteIndex].pixel_count;
			if (sum <= 0)
				continue;

			auto alpha = rgb_table3[paletteIndex].aa;
			pPalette->Entries[paletteIndex] = Color::MakeARGB(alpha, rgb_table3[paletteIndex].rr, rgb_table3[paletteIndex].gg, rgb_table3[paletteIndex].bb);
		}
	}

	bool Dl3Quantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		auto rgb_table3 = make_unique<CUBE3[]>(65536);

		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);
		
		UINT tot_colors = build_table3(rgb_table3.get(), pixels);
		int sqr_tbl[BYTE_MAX + BYTE_MAX + 1];

		for (int i = (-BYTE_MAX); i <= BYTE_MAX; i++)
			sqr_tbl[i + BYTE_MAX] = i * i;

		int* squares3 = &sqr_tbl[BYTE_MAX];

		reduce_table3(rgb_table3.get(), squares3, tot_colors, nMaxColors);
	
		auto qPixels = make_unique<short[]>(bitmapWidth * bitmapHeight);
		quantize_image3(pixels.data(), rgb_table3.get(), nMaxColors, qPixels.get(), squares3, bitmapWidth, bitmapHeight, dither);
		closestMap.clear();
		rightMatches.clear();

		auto pPaletteBytes = make_unique<byte[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		GetQuantizedPalette(pPalette, rgb_table3.get());
		rgb_table3.reset();
		
		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if(nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				std::swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}
