#pragma once
/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2021 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnLABQuantizer.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include "BlueNoise.h"
#include "GilbertCurve.h"
#include <ctime>
#include <unordered_map>

namespace PnnLABQuant
{
	double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;
	BYTE alphaThreshold = 0xF;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	double ratio = 1.0;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
	unordered_map<ARGB, vector<unsigned short> > closestMap;
	unordered_map<ARGB, unsigned short> nearestMap;

	static const float coeffs[3][3] = {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};

	struct pnnbin {
		float ac = 0, Lc = 0, Ac = 0, Bc = 0, err = 0;
		float cnt = 0;
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

	void find_nn(pnnbin* bins, int idx, bool texicab)
	{
		int nn = 0;
		double err = 1e100;

		auto& bin1 = bins[idx];
		auto n1 = bin1.cnt;
		CIELABConvertor::Lab lab1;
		lab1.alpha = bin1.ac, lab1.L = bin1.Lc, lab1.A = bin1.Ac, lab1.B = bin1.Bc;
		for (int i = bin1.fw; i; i = bins[i].fw) {
			auto n2 = bins[i].cnt;
			auto nerr2 = (n1 * n2) / (n1 + n2);
			if (nerr2 >= err)
				continue;

			CIELABConvertor::Lab lab2;
			lab2.alpha = bins[i].ac, lab2.L = bins[i].Lc, lab2.A = bins[i].Ac, lab2.B = bins[i].Bc;
			auto alphaDiff = hasSemiTransparency ? sqr(lab2.alpha - lab1.alpha) / exp(1.5) : 0;
			auto nerr = nerr2 * alphaDiff;
			if (nerr >= err)
				continue;

			if (hasSemiTransparency || !texicab) {
				nerr += (1 - ratio) * nerr2 * sqr(lab2.L - lab1.L);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * sqr(lab2.A - lab1.A);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * sqr(lab2.B - lab1.B);
			}
			else {
				nerr += (1 - ratio) * nerr2 * abs(lab2.L - lab1.L);
				if (nerr >= err)
					continue;

				nerr += (1 - ratio) * nerr2 * _sqrt(sqr(lab2.A - lab1.A) + sqr(lab2.B - lab1.B));
			}

			if (nerr >= err)
				continue;

			auto deltaL_prime_div_k_L_S_L = CIELABConvertor::L_prime_div_k_L_S_L(lab1, lab2);
			nerr += ratio * nerr2 * sqr(deltaL_prime_div_k_L_S_L);
			if (nerr >= err)
				continue;

			double a1Prime, a2Prime, CPrime1, CPrime2;
			auto deltaC_prime_div_k_L_S_L = CIELABConvertor::C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
			nerr += ratio * nerr2 * sqr(deltaC_prime_div_k_L_S_L);
			if (nerr >= err)
				continue;

			double barCPrime, barhPrime;
			auto deltaH_prime_div_k_L_S_L = CIELABConvertor::H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
			nerr += ratio * nerr2 * sqr(deltaH_prime_div_k_L_S_L);
			if (nerr >= err)
				continue;

			nerr += ratio * nerr2 * CIELABConvertor::R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			if (nerr >= err)
				continue;

			err = nerr;
			nn = i;
		}
		bin1.err = err;
		bin1.nn = nn;
	}

	typedef float (*QuanFn)(const float& cnt, const bool isBlack);
	QuanFn getQuanFn(const UINT& nMaxColors, const short quan_rt) {
		if (quan_rt > 0) {
			if (quan_rt > 1)
				return[](const float& cnt, const bool isBlack) { return (float)(int)pow(cnt, 0.75); };
			if (nMaxColors < 64)
				return[](const float& cnt, const bool isBlack) {
					if(isBlack)
						return (float)(int)pow(cnt, 0.75);
					return (float)(int)_sqrt(cnt);
				};
			return[](const float& cnt, const bool isBlack) {
				if (isBlack)
					return (float)pow(cnt, 0.75);
				return (float)_sqrt(cnt);
			};
		}
		return[](const float& cnt, const bool isBlack) { return cnt; };
	}

	void PnnLABQuantizer::pnnquan(const vector<ARGB>& pixels, ColorPalette* pPalette, UINT& nMaxColors)
	{
		short quan_rt = 1;
		vector<pnnbin> bins(USHRT_MAX + 1);

		/* Build histogram */
		for (int i = 0; i < pixels.size(); ++i) {
			const auto& pixel = pixels[i];
			Color c(pixel);
			if (c.GetA() <= alphaThreshold)
				c = m_transparentColor;

			int index = GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);

			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
			auto& tb = bins[index];
			tb.ac += c.GetA();
			tb.Lc += lab1.L;
			tb.Ac += lab1.A;
			tb.Bc += lab1.B;
			tb.cnt += 1.0;
		}

		/* Cluster nonempty bins at one end of array */
		int maxbins = 0;

		for (int i = 0; i < bins.size(); ++i) {
			if (bins[i].cnt <= 0.0)
				continue;

			auto d = 1.0f / bins[i].cnt;
			bins[i].ac *= d;
			bins[i].Lc *= d;
			bins[i].Ac *= d;
			bins[i].Bc *= d;

			bins[maxbins++] = bins[i];
		}

		auto proportional = sqr(nMaxColors) / maxbins;
		if ((m_transparentPixelIndex >= 0 || hasSemiTransparency) && nMaxColors < 32)
			quan_rt = -1;

		auto weight = min(0.9, nMaxColors * 1.0 / maxbins);
		if (weight > .0015 && weight < .002)
			quan_rt = 2;
		if (weight < .025 && PG < 1) {
			auto delta = 3 * (.025 + weight);
			PG -= delta;
			PB += delta;
			if (nMaxColors >= 64)
				quan_rt = 0;
		}

		if (pixelMap.size() <= nMaxColors) {
			/* Fill palette */
			nMaxColors = pPalette->Count = pixelMap.size();
			int k = 0;
			for (const auto& [pixel, lab] : pixelMap) {
				pPalette->Entries[k] = pixel;

				Color c(pPalette->Entries[k]);
				if (k > 0 && c.GetA() == 0)
					swap(pPalette->Entries[k], pPalette->Entries[0]);
				++k;
			}

			return;
		}

		auto quanFn = getQuanFn(nMaxColors, quan_rt);

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;

			bins[j].cnt = quanFn(bins[j].cnt, j == 0);
		}
		bins[j].cnt = quanFn(bins[j].cnt, j == 0);

		const bool texicab = proportional > .025;		
		
		if (quan_rt != 0 && nMaxColors < 64) {
			if (proportional > .018 && proportional < .022)
				ratio = min(1.0, proportional + weight * exp(3.872));
			else if (proportional > .1)
				ratio = min(1.0, 1.0 - weight);
			else if (proportional > .04)
				ratio = min(1.0, weight * exp(2.28));
			else if (proportional > .03)
				ratio = min(1.0, weight * exp(3.275));
			else
				ratio = min(1.0, proportional - weight * exp(1.997));
		}
		else if (nMaxColors > 256)
			ratio = min(1.0, 1 - 1.0 / proportional);
		else
			ratio = min(1.0, max(.98, 1 - weight * .7));

		if (quan_rt < 0)
			ratio = min(1.0, weight * exp(1.997));

		int h, l, l2;
		/* Initialize nearest neighbors and build heap of them */
		auto heap = make_unique<int[]>(bins.size() + 1);
		for (int i = 0; i < maxbins; ++i) {
			find_nn(bins.data(), i, texicab);
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

		if (quan_rt > 0 && nMaxColors < 64 && (proportional < .023 || proportional > .05) && proportional < .1)
			ratio = min(1.0, proportional - weight * exp(2.347));

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
					find_nn(bins.data(), b1, texicab && proportional < 1);
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
			auto d = 1.0 / (n1 + n2);
			tb.ac = d * (n1 * tb.ac + n2 * nb.ac);
			tb.Lc = d * (n1 * tb.Lc + n2 * nb.Lc);
			tb.Ac = d * (n1 * tb.Ac + n2 * nb.Ac);
			tb.Bc = d * (n1 * tb.Bc + n2 * nb.Bc);
			tb.cnt += nb.cnt;
			tb.mtm = ++i;

			/* Unchain deleted bin */
			bins[nb.bk].fw = nb.fw;
			bins[nb.fw].bk = nb.bk;
			nb.mtm = USHRT_MAX;
		}

		/* Fill palette */
		short k = 0;
		for (int i = 0;; ++k) {
			CIELABConvertor::Lab lab1;
			lab1.alpha = (hasSemiTransparency || m_transparentPixelIndex > -1) ? rint(bins[i].ac) : BYTE_MAX;
			lab1.L = bins[i].Lc, lab1.A = bins[i].Ac, lab1.B = bins[i].Bc;
			pPalette->Entries[k] = CIELABConvertor::LAB2RGB(lab1);
			if (m_transparentPixelIndex >= 0 && lab1.alpha == 0) {
				swap(pPalette->Entries[0], pPalette->Entries[k]);
				pPalette->Entries[0] = m_transparentColor;
			}

			if (!(i = bins[i].fw))
				break;
		}

		if (k < nMaxColors - 1)
			pPalette->Count = nMaxColors = k + 1;
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			c = m_transparentColor;

		double mindist = INT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		const auto nMaxColors = pPalette->Count;
		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette->Entries[i]);
			auto curdist = hasSemiTransparency ? sqr(c2.GetA() - c.GetA()) / exp(1.5) : 0;
			if (curdist > mindist)
				continue;

			getLab(c2, lab2);
			if (nMaxColors <= 4) {
				curdist += sqr(c2.GetR() - c.GetR());
				if (curdist > mindist)
					continue;

				curdist += sqr(c2.GetG() - c.GetG());
				if (curdist > mindist)
					continue;

				curdist += sqr(c2.GetB() - c.GetB());
				if (hasSemiTransparency) {
					if (curdist > mindist)
						continue;
					curdist += sqr(c2.GetA() - c.GetA());
				}
			}
			else if (hasSemiTransparency) {
				curdist += sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.B - lab1.B);
			}
			else if (nMaxColors > 32) {
				curdist += abs(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += _sqrt(sqr(lab2.A - lab1.A) + sqr(lab2.B - lab1.B));
			}
			else {
				auto deltaL_prime_div_k_L_S_L = CIELABConvertor::L_prime_div_k_L_S_L(lab1, lab2);
				curdist += sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double a1Prime, a2Prime, CPrime1, CPrime2;
				auto deltaC_prime_div_k_L_S_L = CIELABConvertor::C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double barCPrime, barhPrime;
				auto deltaH_prime_div_k_L_S_L = CIELABConvertor::H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist += sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				curdist += CIELABConvertor::R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
			}

			if (curdist > mindist)
				continue;
			mindist = curdist;
			k = i;
		}
		nearestMap[argb] = k;
		return k;
	}

	unsigned short closestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos)
	{
		UINT k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return nearestColorIndex(pPalette, argb, pos);

		const auto nMaxColors = pPalette->Count;
		vector<unsigned short> closest(4);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = USHRT_MAX;
			
			int start = 0;
			if(BlueNoise.RAW_BLUE_NOISE[pos & 63] > 6)
				start = 1;
			
			for (; k < nMaxColors; ++k) {
				Color c2(pPalette->Entries[k]);				
				
				auto err = PR * (1 - ratio) * sqr(c2.GetR() - c.GetR());
				if (err >= closest[3])
					continue;

				err += PG * (1 - ratio) * sqr(c2.GetG() - c.GetG());
				if (err >= closest[3])
					continue;

				err += PB * (1 - ratio) * sqr(c2.GetB() - c.GetB());
				if (err >= closest[3])
					continue;

				if (hasSemiTransparency)
					err += PA * (1 - ratio) * sqr(c2.GetA() - c.GetA());
				else {
					for (int i = start; i < 3; ++i) {
						err += ratio * sqr(coeffs[i][0] * (c2.GetR() - c.GetR()));
						if (err >= closest[3])
							break;
						
						err += ratio * sqr(coeffs[i][1] * (c2.GetG() - c.GetG()));
						if (err >= closest[3])
							break;
						
						err += ratio * sqr(coeffs[i][2] * (c2.GetB() - c.GetB()));
						if (err >= closest[3])
							break;
					}
				}

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

		auto MAX_ERR = pPalette->Count;
		if(hasSemiTransparency && MAX_ERR > 32) {
			MAX_ERR <<= 1;
			if (c.GetR() > 0xF0 && c.GetG() > 0xF0 && c.GetB() > 0xF0)
				MAX_ERR >>= 1;
		}

		int idx = 1;
		if (closest[2] == 0 || (rand() % (int)ceil(closest[3] + closest[2])) <= closest[3])
			idx = 0;

		if (closest[idx + 2] >= MAX_ERR)
			return nearestColorIndex(pPalette, argb, pos);
		return closest[idx];
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 64) ? nearestColorIndex : closestColorIndex;
		if (dither)
			return dither_image(pixels, pPalette, ditherFn, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i)
				qPixels[pixelIndex++] = ditherFn(pPalette, pixels[pixelIndex], i + j);
		}
		return true;
	}

	bool PnnLABQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapWidth = pSource->GetWidth();
		const auto bitmapHeight = pSource->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		int semiTransCount = 0;
		GrabPixels(pSource, pixels, semiTransCount, m_transparentPixelIndex, m_transparentColor, alphaThreshold, nMaxColors);
		hasSemiTransparency = semiTransCount > 0;

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

		if (nMaxColors <= 32)
			PR = PG = PB = PA = 1;

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		if ((semiTransCount * 1.0 / pixels.size()) > .099)
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, closestColorIndex, GetColorIndex, qPixels.get(), 1.5f);
		else if (nMaxColors <= 32)
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, closestColorIndex, GetColorIndex, qPixels.get(), 1.25f);
		else {
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, closestColorIndex, GetColorIndex, qPixels.get());
			if (nMaxColors > 256) {
				auto qHPixels = make_unique<ARGB[]>(pixels.size());
				for (int i = 0; i < pixels.size(); ++i) {
					Color c(pPalette->Entries[qPixels[i]]);
					qHPixels[i] = hasSemiTransparency ? c.GetValue() : GetARGBIndex(c, false, m_transparentPixelIndex >= 0);
				}

				pixelMap.clear();
				closestMap.clear();
				nearestMap.clear();
				return ProcessImagePixels(pDest, qHPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
			}
		}

		if (!dither) {
			const auto delta = sqr(nMaxColors) / pixelMap.size();
			auto weight = delta > 0.023 ? 1.0f : (float)(36.921 * delta + 0.906);
			BlueNoise::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, closestColorIndex, GetColorIndex, qPixels.get(), weight);
		}

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		pixelMap.clear();
		closestMap.clear();
		nearestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get(), m_transparentPixelIndex >= 0);
	}

}
