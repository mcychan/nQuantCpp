/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2024 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnQuantizer.h"
#include "GilbertCurve.h"
#include "BlueNoise.h"
#include <unordered_map>

namespace PnnQuant
{
	BYTE alphaThreshold = 0xF;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	double ratio = .5, weight = 1.0;
	ARGB m_transparentColor = Color::Transparent;
	double PR = .299, PG = .587, PB = .114, PA = .3333;
	unordered_map<ARGB, vector<unsigned short> > closestMap;
	unordered_map<ARGB, unsigned short > nearestMap;

	static const float coeffs[3][3] = {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};

	struct pnnbin {
		float ac = 0, rc = 0, gc = 0, bc = 0, err = 0;
		float cnt = 0;
		int nn = 0, fw = 0, bk = 0, tm = 0, mtm = 0;
	};

	void find_nn(pnnbin* bins, int idx)
	{
		int nn = 0;
		float err = 1e100;

		auto& bin1 = bins[idx];
		auto n1 = bin1.cnt;
		auto wa = bin1.ac;
		auto wr = bin1.rc;
		auto wg = bin1.gc;
		auto wb = bin1.bc;

		int start = 0;
		if (BlueNoise::TELL_BLUE_NOISE[idx & 4095] > -88)
			start = (PG < coeffs[0][1]) ? 3 : 1;

		for (int i = bin1.fw; i; i = bins[i].fw) {
			auto n2 = bins[i].cnt, nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;

			auto nerr = 0.0;
			if (hasSemiTransparency) {
				start = 1;
				nerr += nerr2 * (1 - ratio) * PA * sqr(bins[i].ac - wa);
				if (nerr >= err)
					continue;
			}

			nerr += nerr2 * (1 - ratio) * PR * sqr(bins[i].rc - wr);
			if (nerr >= err)
				continue;

			nerr += nerr2 * (1 - ratio) * PG * sqr(bins[i].gc - wg);
			if (nerr >= err)
				continue;

			nerr += nerr2 * (1 - ratio) * PB* sqr(bins[i].bc - wb);
			if (nerr >= err)
				continue;

			for (int j = start; j < 3; ++j) {
				nerr += nerr2 * ratio * sqr(coeffs[j][0] * (bins[i].rc - wr));
				if (nerr >= err)
					break;

				nerr += nerr2 * ratio * sqr(coeffs[j][1] * (bins[i].gc - wg));
				if (nerr >= err)
					break;

				nerr += nerr2 * ratio * sqr(coeffs[j][2] * (bins[i].bc - wb));
				if (nerr >= err)
					break;
			}

			if (nerr >= err)
				continue;
			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	typedef float (*QuanFn)(const float& cnt);
	QuanFn getQuanFn(const UINT& nMaxColors, const short quan_rt) {
		if (quan_rt > 0) {
			if (nMaxColors < 64)
				return[](const float& cnt) {
					return (float)(int)_sqrt(cnt);
				};
			return[](const float& cnt) {
				return (float)_sqrt(cnt);
			};
		}
		if (quan_rt < 0)
			return[](const float& cnt) { return (float)(int)cbrt(cnt); };
		return[](const float& cnt) { return cnt; };
	}

	void pnnquan(const vector<ARGB>& pixels, ARGB* pPalette, UINT& nMaxColors)
	{
		short quan_rt = 1;
		vector<pnnbin> bins(USHRT_MAX + 1);

		/* Build histogram */
		for (const auto& pixel : pixels) {
			Color c(pixel);
			if (c.GetA() <= alphaThreshold)
				c = m_transparentColor;

			int index = GetARGBIndex(c, hasSemiTransparency, nMaxColors < 64 || m_transparentPixelIndex >= 0);
			auto& tb = bins[index];
			tb.ac += c.GetA();
			tb.rc += c.GetR();
			tb.gc += c.GetG();
			tb.bc += c.GetB();
			tb.cnt += 1.0;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < bins.size(); ++i) {
			if (bins[i].cnt <= 0)
				continue;

			auto& tb = bins[i];
			double d = 1.0 / tb.cnt;
			tb.ac *= d;
			tb.rc *= d;
			tb.gc *= d;
			tb.bc *= d;
			
			bins[maxbins++] = tb;
		}

		if (nMaxColors < 16)
			quan_rt = -1;

		weight = min(0.9, nMaxColors * 1.0 / maxbins);
		if (weight < .03 && PG < 1 && PG >= coeffs[0][1]) {
			PR = PG = PB = PA = 1;
			if (nMaxColors >= 64)
				quan_rt = 0;
		}

		auto quanFn = getQuanFn(nMaxColors, quan_rt);

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;

			bins[j].cnt = quanFn(bins[j].cnt);
		}
		bins[j].cnt = quanFn(bins[j].cnt);

		auto heap = make_unique<int[]>(bins.size() + 1);
		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins.data(), i);
			/* Push slot on heap */
			auto err = bins[i].err;
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
				if (tb.mtm == USHRT_MAX) /* Deleted node */
					b1 = heap[1] = heap[heap[0]--];
				else /* Too old error value */
				{
					find_nn(bins.data(), b1);
					tb.tm = i;
				}
				/* Push slot down */
				auto err = bins[b1].err;
				for (l = 1; (l2 = l + l) <= heap[0]; l = l2) {
					if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
						++l2;
					if (err <= bins[h = heap[l2]].err)
						break;
					heap[l] = h;
				}
				heap[l] = b1;
			}

			/* Do a merge */
			auto& tb = bins[b1];
			auto& nb = bins[tb.nn];
			auto n1 = tb.cnt;
			auto n2 = nb.cnt;
			auto d = 1.0f / (n1 + n2);
			tb.ac = d * rint(n1 * tb.ac + n2 * nb.ac);
			tb.rc = d * rint(n1 * tb.rc + n2 * nb.rc);
			tb.gc = d * rint(n1 * tb.gc + n2 * nb.gc);
			tb.bc = d * rint(n1 * tb.bc + n2 * nb.bc);
			tb.cnt += n2;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = USHRT_MAX;
		}

		/* Fill palette */
		UINT k = 0;
		for (int i = 0;; ++k) {
			auto alpha = (hasSemiTransparency || m_transparentPixelIndex > -1) ? rint(bins[i].ac) : BYTE_MAX;
			pPalette[k] = Color::MakeARGB(alpha, (int) bins[i].rc, (int) bins[i].gc, (int) bins[i].bc);

			if (!(i = bins[i].fw))
				break;
		}

		if (k < nMaxColors - 1)
			nMaxColors = k + 1;
	}

	unsigned short nearestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			c = m_transparentColor;

		if (nMaxColors > 2 && m_transparentPixelIndex >= 0 && c.GetA() > alphaThreshold)
			k = 1;
		
		auto pr = PR, pg = PG, pb = PB, pa = PA;
		if(nMaxColors < 3)
			pr = pg = pb = pa = 1;

		double mindist = INT_MAX;		
		for (UINT i = k; i < nMaxColors; ++i) {
			Color c2(pPalette[i]);
			double curdist = pa * sqr(c2.GetA() - c.GetA());
			if (curdist > mindist)
				continue;

			curdist += pr * sqr(c2.GetR() - c.GetR());
			if (curdist > mindist)
				continue;

			curdist += pg * sqr(c2.GetG() - c.GetG());
			if (curdist > mindist)
				continue;

			curdist += pb * sqr(c2.GetB() - c.GetB());
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		nearestMap[argb] = k;
		return k;
	}

	unsigned short closestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		UINT k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);

		vector<unsigned short> closest(4);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = USHRT_MAX;

			auto pr = PR, pg = PG, pb = PB, pa = PA;
			if(nMaxColors < 3)
				pr = pg = pb = pa = 1;

			for (; k < nMaxColors; ++k) {
				Color c2(pPalette[k]);
				auto err = pr * sqr(c2.GetR() - c.GetR());
				if (err >= closest[3])
					break;

				err += pg* sqr(c2.GetG() - c.GetG());
				if (err >= closest[3])
					break;

				err += pb* sqr(c2.GetB() - c.GetB());
				if (err >= closest[3])
					break;

				if (hasSemiTransparency)
					err += pa * sqr(c2.GetA() - c.GetA());

				if (err < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = err;
				}
				else if (err < closest[3]) {
					closest[1] = k;
					closest[3] = err;
				}
			}

			if (closest[3] == USHRT_MAX)
				closest[1] = closest[0];

			closestMap[argb] = closest;
		}
		else
			closest = got->second;

		auto MAX_ERR = nMaxColors << 2;
		int idx = (pos + 1) % 2;
		if (closest[3] * .67 < (closest[3] - closest[2]))
			idx = 0;
		else if (closest[0] > closest[1])
			idx = pos % 2;

		if (closest[idx + 2] >= MAX_ERR || (m_transparentPixelIndex >= 0 && closest[idx] == 0))
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);
		return closest[idx];
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{		
		if (dither) 
			return dither_image(pixels, pPalette->Entries, nMaxColors, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i, ++pixelIndex)
				qPixels[pixelIndex] = ditherFn(pPalette->Entries, nMaxColors, pixels[pixelIndex], i + j);
		}

		BlueNoise::dither(width, height, pixels, pPalette->Entries, nMaxColors, ditherFn, GetColorIndex, qPixels);
		return true;
	}

	bool PnnQuantizer::QuantizeImage(const vector<ARGB>& pixels, const UINT bitmapWidth, ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		if (nMaxColors <= 32)
			PR = PG = PB = PA = 1;
		else {
			PR = coeffs[0][0]; PG = coeffs[0][1]; PB = coeffs[0][2];
		}

		const auto bitmapHeight = pixels.size() / bitmapWidth;

		if (nMaxColors > 2)
			pnnquan(pixels, pPalette, nMaxColors);
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette[0] = m_transparentColor;
				pPalette[1] = Color::Black;
			}
			else {
				pPalette[0] = Color::Black;
				pPalette[1] = Color::White;
			}
		}

		DitherFn ditherFn = dither ? nearestColorIndex : closestColorIndex;
		if (hasSemiTransparency)
			weight *= -1;

		vector<float> saliencies;

		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get(), saliencies.data(), weight);

			closestMap.clear();
			nearestMap.clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get(), saliencies.data(), weight);

		if (!dither)
			BlueNoise::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get());

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette[k] = m_transparentColor;
			else if (pPalette[k] != m_transparentColor)
				swap(pPalette[0], pPalette[1]);
		}
		closestMap.clear();
		nearestMap.clear();

		return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex >= 0);
	}

	bool PnnQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapWidth = pSource->GetWidth();
		const auto bitmapHeight = pSource->GetHeight();
		const auto area = (size_t) (bitmapWidth * bitmapHeight);

		vector<ARGB> pixels(area);
		int semiTransCount = 0;
		GrabPixels(pSource, pixels, semiTransCount, m_transparentPixelIndex, m_transparentColor, alphaThreshold, nMaxColors);
		hasSemiTransparency = semiTransCount > 0;
		
		if (nMaxColors > 256) {
			auto pPalettes = make_unique<ARGB[]>(nMaxColors);
			auto pPalette = pPalettes.get();
			return QuantizeImage(pixels, bitmapWidth, pPalette, pDest, nMaxColors, dither);
		}

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		auto result = QuantizeImage(pixels, bitmapWidth, pPalette->Entries, pDest, nMaxColors, dither);
		pDest->SetPalette(pPalette);
		return result;
	}

}
