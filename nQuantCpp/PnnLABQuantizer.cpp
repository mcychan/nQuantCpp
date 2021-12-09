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
	double PR = .2126, PG = .7152, PB = .0722;
	BYTE alphaThreshold = 0;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	double ratio = 1.0;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
	unordered_map<ARGB, vector<unsigned short> > closestMap;
	unordered_map<ARGB, unsigned short> nearestMap;

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

	void find_nn(pnnbin* bins, int idx)
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
			auto alphaDiff = hasSemiTransparency ? abs(lab2.alpha - lab1.alpha) : 0;
			auto nerr = nerr2 * sqr(alphaDiff) / exp(1.5);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * sqr(lab2.L - lab1.L);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * sqr(lab2.A - lab1.A);
			if (nerr >= err)
				continue;

			nerr += (1 - ratio) * nerr2 * sqr(lab2.B - lab1.B);

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

	typedef float (*QuanFn)(const float& cnt);
	QuanFn getQuanFn(const UINT& nMaxColors, const short quan_rt) {
		if (quan_rt > 0) {
			if (quan_rt > 1)
				return[](const float& cnt) { return (float)(int) pow(cnt, 0.75); };
			if (nMaxColors < 64)
				return[](const float& cnt) { return (float)(int) _sqrt(cnt); };
			return[](const float& cnt) { return (float) _sqrt(cnt); };
		}
		return[](const float& cnt) { return cnt; };
	}

	void PnnLABQuantizer::pnnquan(const vector<ARGB>& pixels, ColorPalette* pPalette, UINT& nMaxColors, short quan_rt)
	{
		vector<pnnbin> bins(USHRT_MAX + 1);

		/* Build histogram */
		for (const auto& pixel : pixels) {
			// !!! Can throw gamma correction in here, but what to do about perceptual
			// !!! nonuniformity then?			
			Color c(pixel);

			int index = GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);

			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
			bins[index].ac += c.GetA();
			bins[index].Lc += lab1.L;
			bins[index].Ac += lab1.A;
			bins[index].Bc += lab1.B;
			bins[index].cnt += 1.0;
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
		
		auto weight = nMaxColors * 1.0 / maxbins;
		if (weight > .0015 && weight < .002)
			quan_rt = 2;

		auto quanFn = getQuanFn(nMaxColors, quan_rt);

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;

			bins[j].cnt = quanFn(bins[j].cnt);			
		}
		bins[j].cnt = quanFn(bins[j].cnt);
		
		int h, l, l2;
		if (quan_rt != 0 && nMaxColors < 64) {
			if (proportional > .018 && proportional < .022)
				ratio = min(1.0, proportional + nMaxColors * exp(3.872) / maxbins);
			else if (proportional > .1)
				ratio = min(1.0, proportional - nMaxColors * exp(3.23) / maxbins);
			else
				ratio = min(1.0, proportional - nMaxColors * exp(1.997) / maxbins);
		}
		else if (nMaxColors > 256)
			ratio = min(hasSemiTransparency ? 0.0 : 1.0, 1 - 1.0 / proportional);
		else
			ratio = min(hasSemiTransparency ? 0.0 : 1.0, 0.14 * exp(4.681 * proportional));

		if (quan_rt < 0)
			ratio = min(1.0, nMaxColors * exp(1.997) / maxbins);
				
		/* Initialize nearest neighbors and build heap of them */
		auto heap = make_unique<int[]>(bins.size() + 1);
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

		if (quan_rt > 0 && nMaxColors < 64 && (proportional < .023 || proportional > .05))
			ratio = min(1.0, proportional - nMaxColors * exp(2.347) / maxbins);

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
			lab1.alpha = rint(bins[i].ac);
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

	unsigned short nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return k;

		double mindist = USHRT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette->Entries[i]);
			double curdist = hasSemiTransparency ? sqr(c2.GetA() - c.GetA()) : 0;
			if (curdist > mindist)
				continue;
			
			if (nMaxColors > 32 || nMaxColors <= 4 || hasSemiTransparency) {
				curdist += PR * sqr(c2.GetR() - c.GetR());
				if (curdist > mindist)
					continue;

				curdist += PG * sqr(c2.GetG() - c.GetG());
				if (curdist > mindist)
					continue;

				curdist += PB * sqr(c2.GetB() - c.GetB());
				if (PB < 1) {
					if (curdist > mindist)
						continue;

					getLab(c2, lab2);
					curdist += sqr(lab2.B - lab1.B) / 2.0;
				}
			}
			else {
				getLab(c2, lab2);
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

	unsigned short closestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		UINT k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return k;

		vector<unsigned short> closest(4);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = USHRT_MAX;

			for (; k < nMaxColors; ++k) {
				Color c2(pPalette->Entries[k]);		
				auto err = PR * sqr(c2.GetR() - c.GetR()) + PG * sqr(c2.GetG() - c.GetG()) + PB * sqr(c2.GetB() - c.GetB());
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
		if (closest[2] == 0 || (rand() % (int)ceil(closest[3] + closest[2])) <= closest[3]) {
			if (closest[2] > MAX_ERR)
				return nearestColorIndex(pPalette, nMaxColors, argb);
			return closest[0];
		}

		if (closest[3] > MAX_ERR)
			return nearestColorIndex(pPalette, nMaxColors, argb);
		return closest[1];
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

		for (int pixelIndex = 0; pixelIndex < (width * height); ++pixelIndex)
			qPixels[pixelIndex] = ditherFn(pPalette, nMaxColors, pixels[pixelIndex]);		
		return true;
	}	

	bool PnnLABQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapWidth = pSource->GetWidth();
		const auto bitmapHeight = pSource->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors > 2)
			pnnquan(pixels, pPalette, nMaxColors, 1);
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

		bool noBias = m_transparentPixelIndex >= 0 || nMaxColors < 64;
		if (noBias)
			PR = PG = PB = 1;
		else if (bitmapWidth < 512 || bitmapHeight < 512) {
			PR = 0.299; PG = 0.587; PB = 0.114;
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		DitherFn ditherFn = hasSemiTransparency ? nearestColorIndex : closestColorIndex;
		if (hasSemiTransparency && nMaxColors <= 256)
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, ditherFn, GetColorIndex, qPixels.get(), 1.75f);
		else if (nMaxColors < 64 && nMaxColors > 32)
			quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);
		else if(nMaxColors <= 32)
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, ditherFn, GetColorIndex, qPixels.get(), 1.5f);
		else {
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, ditherFn, GetColorIndex, qPixels.get());
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
			BlueNoise::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, ditherFn, GetColorIndex, qPixels.get(), weight);
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
