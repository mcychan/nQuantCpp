/* Fast pairwise nearest neighbor based algorithm for multilevel thresholding
Copyright (C) 2004-2016 Mark Tyler and Dmitry Groshev
Copyright (c) 2018-2024 Miller Cy Chan
* error measure; time used is proportional to number of bins squared - WJ */

#include "stdafx.h"
#include "PnnLABQuantizer.h"
#include "bitmapUtilities.h"
#include "BlueNoise.h"
#include "GilbertCurve.h"
#include <ctime>

namespace PnnLABQuant
{
	double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;
	BYTE alphaThreshold = 0xF;
	double weight = 1.0;

	ARGB m_transparentColor = Color::Transparent;

	static const float coeffs[3][3] = {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};

	PnnLABQuantizer::PnnLABQuantizer() {
	}

	PnnLABQuantizer::PnnLABQuantizer(const PnnLABQuantizer& quantizer) {
		hasSemiTransparency = quantizer.hasSemiTransparency;
		m_transparentPixelIndex = quantizer.m_transparentPixelIndex;
		saliencies = quantizer.saliencies;
		pixelMap.insert(quantizer.pixelMap.begin(), quantizer.pixelMap.end());
		isGA = true;
		proportional = quantizer.proportional;
	}

	void PnnLABQuantizer::getLab(const Color& c, CIELABConvertor::Lab& lab1)
	{
		auto got = pixelMap.find(c.GetValue());
		if (got == pixelMap.end()) {
			CIELABConvertor::RGB2LAB(c, lab1);
			pixelMap[c.GetValue()] = lab1;
		}
		else
			lab1 = got->second;
	}

	void PnnLABQuantizer::find_nn(pnnbin* bins, int idx, bool texicab)
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

	typedef float (*QuanFn)(const float& cnt);
	QuanFn getQuanFn(const UINT& nMaxColors, const short quan_rt) {
		if (quan_rt > 0) {
			if (quan_rt > 1)
				return[](const float& cnt) { return (float)(int)pow(cnt, 0.75); };
			if (nMaxColors < 64)
				return[](const float& cnt) {
					return (float)(int)_sqrt(cnt);
				};
			return[](const float& cnt) {
				return (float)_sqrt(cnt);
			};
		}
		return[](const float& cnt) { return cnt; };
	}

	void PnnLABQuantizer::pnnquan(const vector<ARGB>& pixels, ARGB* pPalette, UINT& nMaxColors)
	{
		short quan_rt = 1;
		vector<pnnbin> bins(USHRT_MAX + 1);
		saliencies.resize(pixels.size());
		auto saliencyBase = .1f;

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
			if(lab1.alpha > alphaThreshold && nMaxColors < 32)
				saliencies[i] = saliencyBase + (1 - saliencyBase) * lab1.L / 100.0f;
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

		proportional = sqr(nMaxColors) / maxbins;
		if ((m_transparentPixelIndex >= 0 || hasSemiTransparency) && nMaxColors < 32)
			quan_rt = -1;

		weight = min(0.9, nMaxColors * 1.0 / maxbins);
		if (weight > .0015 && weight < .002)
			quan_rt = 2;
		if (weight < .04 && PG < 1 && PG >= coeffs[0][1]) {
			auto delta = exp(1.75) * weight;
			PG -= delta;
			PB += delta;
			if (nMaxColors >= 64)
				quan_rt = 0;
		}

		if (pixelMap.size() <= nMaxColors) {
			/* Fill palette */
			nMaxColors = pixelMap.size();
			int k = 0;
			for (const auto& [pixel, lab] : pixelMap) {
				pPalette[k] = pixel;

				Color c(pPalette[k]);
				if (k > 0 && c.GetA() == 0)
					swap(pPalette[k], pPalette[0]);
				++k;
			}

			return;
		}

		auto quanFn = getQuanFn(nMaxColors, quan_rt);

		int j = 0;
		for (; j < maxbins - 1; ++j) {
			bins[j].fw = j + 1;
			bins[j + 1].bk = j;

			bins[j].cnt = quanFn(bins[j].cnt);
		}
		bins[j].cnt = quanFn(bins[j].cnt);

		const bool texicab = proportional > .025;
		
		if(!isGA) {
			if (hasSemiTransparency)
				ratio = .5;
			else if (quan_rt != 0 && nMaxColors < 64) {
				if (proportional > .018 && proportional < .022)
					ratio = min(1.0, proportional + weight * exp(3.872));
				else if (proportional > .1)
					ratio = min(1.0, 1.0 - weight);
				else if (proportional > .04)
					ratio = min(1.0, weight * exp(2.28));
				else if (proportional > .03)
					ratio = min(1.0, weight * exp(3.275));
				else {
					auto beta = (maxbins % 2 == 0) ? -1 : 1;
					ratio = min(1.0, proportional + beta * weight * exp(1.997));
				}
			}
			else if (nMaxColors > 256)
				ratio = min(1.0, 1 - 1.0 / proportional);
			else
				ratio = min(1.0, max(.98, 1 - weight * .7));
	
			if (!hasSemiTransparency && quan_rt < 0)
				ratio = min(1.0, weight * exp(1.997));
		}

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

		if (!isGA && quan_rt > 0 && nMaxColors < 64 && (proportional < .023 || proportional > .05) && proportional < .1)
			ratio = min(1.0, proportional - weight * exp(2.347));
		else if (isGA)
			ratio = ratioY;

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
			tb.cnt += n2;
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
			pPalette[k] = CIELABConvertor::LAB2RGB(lab1);

			if (!(i = bins[i].fw))
				break;
		}

		if (k < nMaxColors - 1)
			nMaxColors = k + 1;
	}

	unsigned short PnnLABQuantizer::nearestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			c = m_transparentColor;

		if (nMaxColors > 2 && hasAlpha() && c.GetA() > alphaThreshold)
			k = 1;

		double mindist = INT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);
		
		for (UINT i = k; i < nMaxColors; ++i) {
			Color c2(pPalette[i]);
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

	unsigned short PnnLABQuantizer::closestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		UINT k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);

		vector<unsigned short> closest(4);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = USHRT_MAX;
			
			int start = 0;
			if(BlueNoise::TELL_BLUE_NOISE[pos & 4095] > -88)
				start = 1;
			
			for (; k < nMaxColors; ++k) {
				Color c2(pPalette[k]);
				
				auto err = PR * (1 - ratio) * sqr(c2.GetR() - c.GetR());
				if (err >= closest[3])
					continue;

				err += PG * (1 - ratio) * sqr(c2.GetG() - c.GetG());
				if (err >= closest[3])
					continue;

				err += PB * (1 - ratio) * sqr(c2.GetB() - c.GetB());
				if (err >= closest[3])
					continue;

				if (hasSemiTransparency) {
					err += PA * (1 - ratio) * sqr(c2.GetA() - c.GetA());
					start = 1;
				}
				
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

		auto MAX_ERR = nMaxColors;
		if(PG < coeffs[0][1] && BlueNoise::TELL_BLUE_NOISE[pos & 4095] > -88)
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);

		int idx = 1;
		if (closest[2] == 0 || (rand() % (int)ceil(closest[3] + closest[2])) <= closest[3])
			idx = 0;

		if (closest[idx + 2] >= MAX_ERR || (hasAlpha() && closest[idx] == 0))
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);
		return closest[idx];
	}

	void PnnLABQuantizer::clear()
	{
		saliencies.clear();
		closestMap.clear();
		nearestMap.clear();
	}

	bool PnnLABQuantizer::IsGA() const {
		return isGA;
	}

	bool PnnLABQuantizer::hasAlpha() const {
		return m_transparentPixelIndex >= 0;
	}
	
	void PnnLABQuantizer::setRatio(double ratioX, double ratioY) {
		ratio = min(1.0, ratioX);
		this->ratioY = min(1.0, ratioY);
		clear();
	}

	void PnnLABQuantizer::grabPixels(Bitmap* srcImg, vector<ARGB>& pixels, UINT& nMaxColors, bool& hasSemiTransparency)
	{
		int semiTransCount = 0;
		GrabPixels(srcImg, pixels, semiTransCount, m_transparentPixelIndex, m_transparentColor, alphaThreshold, nMaxColors);
		this->hasSemiTransparency = hasSemiTransparency = semiTransCount > 0;
	}
	
	bool PnnLABQuantizer::quantize_image(const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		DitherFn ditherFn;
		if (hasAlpha() || nMaxColors < 64)
			ditherFn = [&](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
				return nearestColorIndex(pPalette, nMaxColors, argb, pos);
		};
		else
			ditherFn = [&](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
			return closestColorIndex(pPalette, nMaxColors, argb, pos);
		};
		
		if (dither)
			return dither_image(pixels, pPalette, nMaxColors, ditherFn, hasSemiTransparency, m_transparentPixelIndex, qPixels, width, height);

		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i)
				qPixels[pixelIndex++] = ditherFn(pPalette, nMaxColors, pixels[pixelIndex], i + j);
		}
		return true;
	}

	bool PnnLABQuantizer::QuantizeImage(const vector<ARGB>& pixels, const UINT bitmapWidth, ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither)
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
		
		if (hasSemiTransparency)
			weight *= -1;

		auto GetColorIndex = [&](const Color& c) -> int {
			return GetARGBIndex(c, hasSemiTransparency, hasAlpha());
		};
		DitherFn ditherFn;
		if (hasAlpha() || nMaxColors < 64)
			ditherFn = [&](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);
		};
		else
			ditherFn = [&](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
			return closestColorIndex(pPalette, nMaxColors, argb, pos);
		};				

		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get(), saliencies.data(), weight);

			pixelMap.clear();
			clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get(), saliencies.data(), weight);

		if (!dither) {
			const auto delta = sqr(nMaxColors) / pixelMap.size();
			auto weight = delta > 0.023 ? 1.0f : (float)(36.921 * delta + 0.906);
			BlueNoise::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, ditherFn, GetColorIndex, qPixels.get(), weight);
		}

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette[k] = m_transparentColor;
			else if (pPalette[k] != m_transparentColor)
				swap(pPalette[0], pPalette[1]);
		}
		pixelMap.clear();
		clear();

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPals = (ColorPalette*) pPaletteBytes.get();
		pPals->Count = nMaxColors;
		for (UINT k = 0; k < nMaxColors; ++k)
			pPals->Entries[k] = pPalette[k];
		return ProcessImagePixels(pDest, pPals, qPixels.get(), m_transparentPixelIndex >= 0);
	}

	bool PnnLABQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapWidth = pSource->GetWidth();
		const auto bitmapHeight = pSource->GetHeight();
		const auto area = (size_t) (bitmapWidth * bitmapHeight);

		vector<ARGB> pixels(area);
		int semiTransCount = 0;
		grabPixels(pSource, pixels, nMaxColors, hasSemiTransparency);		
		
		auto pPalettes = make_unique<ARGB[]>(nMaxColors);
		auto pPalette = pPalettes.get();

		return QuantizeImage(pixels, bitmapWidth, pPalette, pDest, nMaxColors, dither);
	}

}
