/*
 * DL3 Quantization
 * ================
 *
 * Author: Dennis Lee   E-mail: denlee@ecf.utoronto.ca
 *
 * Copyright (C) 1993-1997 Dennis Lee
 * Copyright (c) 2019-2021 Miller Cy Chan
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
#include "BlueNoise.h"
#include <unordered_map>

namespace Dl3Quant
{
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, vector<unsigned short> > closestMap;

	using namespace std;

	struct CUBE3 {
		int a, r, g, b;
		int aa, rr, gg, bb;
		UINT cc;
		UINT pixel_count = 0;
		double err;
	};

	void setARGB(CUBE3& rec)
	{
		UINT v = rec.pixel_count, v2 = v >> 1;
		rec.aa = (rec.a + v2) / v;
		rec.rr = (rec.r + v2) / v;
		rec.gg = (rec.g + v2) / v;
		rec.bb = (rec.b + v2) / v;
	}

	double calc_err(CUBE3* rgb_table3, const int* squares3, const UINT& c1, const UINT& c2)
	{
		auto P1 = rgb_table3[c1].pixel_count;
		auto P2 = rgb_table3[c2].pixel_count;
		auto P3 = P1 + P2;

		if (P3 == 0)
			return UINT_MAX;

		int A3 = (rgb_table3[c1].a + rgb_table3[c2].a + (P3 >> 1)) / P3;
		int R3 = (rgb_table3[c1].r + rgb_table3[c2].r + (P3 >> 1)) / P3;
		int G3 = (rgb_table3[c1].g + rgb_table3[c2].g + (P3 >> 1)) / P3;
		int B3 = (rgb_table3[c1].b + rgb_table3[c2].b + (P3 >> 1)) / P3;

		int A1 = rgb_table3[c1].aa;
		int R1 = rgb_table3[c1].rr;
		int G1 = rgb_table3[c1].gg;
		int B1 = rgb_table3[c1].bb;

		int A2 = rgb_table3[c2].aa;
		int R2 = rgb_table3[c2].rr;
		int G2 = rgb_table3[c2].gg;
		int B2 = rgb_table3[c2].bb;

		double dist1 = squares3[A3 - A1] + squares3[R3 - R1] + squares3[G3 - G1] + squares3[B3 - B1];
		dist1 *= P1;

		double dist2 = squares3[A2 - A3] + squares3[R2 - R3] + squares3[G2 - G3] + squares3[B2 - B3];
		dist2 *= P2;

		return (dist1 + dist2);
	}

	void build_table3(CUBE3* rgb_table3, ARGB argb)
	{
		Color c(argb);
		int index = GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);

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
		for (int i = 0; i < USHRT_MAX + 1; ++i) {
			if (rgb_table3[i].pixel_count > 0) {
				setARGB(rgb_table3[i]);
				rgb_table3[tot_colors] = rgb_table3[i];
				++tot_colors;
			}
		}
		return tot_colors;
	}

	void recount_next(CUBE3* rgb_table3, const int* squares3, const UINT& tot_colors, const UINT& i)
	{
		UINT c2 = 0;
		double err = UINT_MAX;
		for (UINT j = i + 1; j < tot_colors; ++j) {
			auto cur_err = calc_err(rgb_table3, squares3, i, j);
			if (cur_err < err) {
				err = cur_err;
				c2 = j;
			}
		}
		rgb_table3[i].err = err;
		rgb_table3[i].cc = c2;
	}

	void recount_dist(CUBE3* rgb_table3, const int* squares3, const UINT& tot_colors, const UINT& c1)
	{
		recount_next(rgb_table3, squares3, tot_colors, c1);
		for (int i = 0; i < c1; ++i) {
			if (rgb_table3[i].cc == c1)
				recount_next(rgb_table3, squares3, tot_colors, i);
			else {
				auto cur_err = calc_err(rgb_table3, squares3, i, c1);
				if (cur_err < rgb_table3[i].err) {
					rgb_table3[i].err = cur_err;
					rgb_table3[i].cc = c1;
				}
			}
		}
	}

	void reduce_table3(CUBE3* rgb_table3, const int* squares3, UINT tot_colors, const UINT& num_colors)
	{
		UINT i = 0;
		for (; i < (tot_colors - 1); ++i)
			recount_next(rgb_table3, squares3, tot_colors, i);

		rgb_table3[i].err = UINT_MAX;
		rgb_table3[i].cc = tot_colors;

		UINT c1 = 0, grand_total = tot_colors - num_colors;
		while (tot_colors > num_colors) {
			UINT err = UINT_MAX;
			for (i = 0; i < tot_colors; ++i) {
				if (rgb_table3[i].err < err) {
					err = rgb_table3[i].err;
					c1 = i;
				}
			}
			auto c2 = rgb_table3[c1].cc;
			rgb_table3[c2].a += rgb_table3[c1].a;
			rgb_table3[c2].r += rgb_table3[c1].r;
			rgb_table3[c2].g += rgb_table3[c1].g;
			rgb_table3[c2].b += rgb_table3[c1].b;
			rgb_table3[c2].pixel_count += rgb_table3[c1].pixel_count;
			setARGB(rgb_table3[c2]);

			rgb_table3[c1] = rgb_table3[--tot_colors];
			rgb_table3[tot_colors - 1].err = UINT_MAX;
			rgb_table3[tot_colors - 1].cc = tot_colors;

			for (i = 0; i < c1; ++i) {
				if (rgb_table3[i].cc == tot_colors)
					rgb_table3[i].cc = c1;
			}

			for (i = c1 + 1; i < tot_colors; ++i) {
				if (rgb_table3[i].cc == tot_colors)
					recount_next(rgb_table3, squares3, tot_colors, i);
			}

			recount_dist(rgb_table3, squares3, tot_colors, c1);
			if (c2 != tot_colors)
				recount_dist(rgb_table3, squares3, tot_colors, c2);
		}
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos)
	{
		unsigned short k = 0;
		Color c(argb);

		double mindist = INT_MAX;
		const auto nMaxColors = pPalette->Count;
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

	unsigned short closestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos)
	{
		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= 0xF)
			c = m_transparentColor;

		const auto nMaxColors = pPalette->Count;

		vector<unsigned short> closest(5);
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

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither)
			return dither_image(pixels, pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i, ++pixelIndex)
				qPixels[pixelIndex] = ditherFn(pPalette, pixels[pixelIndex], i + j);
		}

		BlueNoise::dither(width, height, pixels, pPalette, ditherFn, GetColorIndex, qPixels);
		return true;
	}

	void GetQuantizedPalette(ColorPalette* pPalette, const CUBE3* rgb_table3)
	{
		for (UINT k = 0; k < pPalette->Count; ++k) {
			UINT sum = rgb_table3[k].pixel_count;
			if (sum > 0) {
				pPalette->Entries[k] = Color::MakeARGB(rgb_table3[k].aa, rgb_table3[k].rr, rgb_table3[k].gg, rgb_table3[k].bb);

				if (m_transparentPixelIndex >= 0 && pPalette->Entries[k] == m_transparentColor)
					swap(pPalette->Entries[0], pPalette->Entries[k]);
			}
		}
	}

	bool Dl3Quantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor, nMaxColors);

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors > 2) {
			auto rgb_table3 = make_unique<CUBE3[]>(USHRT_MAX + 1);
			UINT tot_colors = build_table3(rgb_table3.get(), pixels);
			int sqr_tbl[BYTE_MAX + BYTE_MAX + 1];

			for (int i = (-BYTE_MAX); i <= BYTE_MAX; ++i)
				sqr_tbl[i + BYTE_MAX] = i * i;

			auto squares3 = &sqr_tbl[BYTE_MAX];

			reduce_table3(rgb_table3.get(), squares3, tot_colors, nMaxColors);

			GetQuantizedPalette(pPalette, rgb_table3.get());
		}
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
			closestMap.clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);
		closestMap.clear();

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		return ProcessImagePixels(pDest, pPalette, qPixels.get(), m_transparentPixelIndex >= 0);
	}

}
