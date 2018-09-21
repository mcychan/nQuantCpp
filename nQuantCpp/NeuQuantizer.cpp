#pragma once

#include "stdafx.h"
#include "NeuQuantizer.h"
#include <map>

namespace NeuralNet
{
	//====================
	// NeuralNet Color quantizing

	/* NeuQuant Neural-Net Quantization Algorithm
	* ------------------------------------------
	*
	* Copyright (c) 1994 Anthony Dekker
	*
	* NEUQUANT Neural-Net quantization algorithm by Anthony Dekker, 1994.
	* See "Kohonen neural networks for optimal colour quantization"
	* in "Network: Computation in Neural Systems" Vol. 5 (1994) pp 351-367.
	* for a discussion of the algorithm.
	* See also  http://www.acm.org/~dekker/NEUQUANT.HTML
	*
	* Any party obtaining a copy of these files from the author, directly or
	* indirectly, is granted, free of charge, a full and unrestricted irrevocable,
	* world-wide, paid up, royalty-free, nonexclusive right and license to deal
	* in this software and documentation files (the "Software"), including without
	* limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
	* and/or sell copies of the Software, and to permit persons who receive
	* copies from any such party to do so, with the only requirement being
	* that this copyright notice remain intact.
	*/

	const int ncycles = 115;			// no. of learning cycles

	const int netsize = 256;		// number of colours used
	const int specials = 3;		// number of reserved colours used
	const int bgColour = specials - 1;	// reserved background colour
	const int cutnetsize = netsize - specials;
	const int maxnetpos = netsize - 1;

	const int initrad = netsize >> 3;   // for 256 cols, radius starts at 32
	const double initradius = initrad * 1.0;
	const int radiusbiasshift = 8;
	const int radiusbias = 1 << radiusbiasshift;
	const int initBiasRadius = initrad * radiusbias;
	const int radiusdec = 30; // factor of 1/30 each cycle

	const int alphabiasshift = 10;			// alpha starts at 1
	const int initalpha = 1 << alphabiasshift; // biased by 10 bits
	const int alpharadbshift = alphabiasshift + radiusbiasshift;
	const double alpharadbias = (double)(1 << alpharadbshift);

	/* defs for freq and bias */
	const int gammashift = 10;                  /* gamma = 1024 */
	const double gamma = (double)(1 << gammashift);
	const int betashift = 10;
	const double beta = (1.0 / (1 << betashift));/* beta = 1/1024 */
	const double betagamma = (double)(1 << (gammashift - betashift));

	struct nq_pixel                        /* ABGRc */
	{
		double al, b, g, h, r;
	};

	nq_pixel network[netsize]; // the network itself

	int netindex[netsize] = { 0 }; // for network lookup - really 256
	double biasvalues[netsize] = { 0 };  // Biasvalues: based on frequency of nearest pixels

	double bias[netsize] = { 0.0 };  // bias and freq arrays for learning
	double freq[netsize] = { 0.0 };
	double radpower[initrad] = { 0.0 };

	// four primes near 500 - assume no image has a length so large
	// that it is divisible by all four primes

	const int primes[] = { 499, 491, 487, 503 };
	double gamma_correction = 1.0;         // 1.0/2.2 usually

	bool hasTransparency = false;
	ARGB m_transparentColor;
	vector<ARGB> pixels;
	map<ARGB, vector<short> > closestMap;

	inline double biasvalue(unsigned int temp)
	{
		return biasvalues[temp];
	}

	inline double colorimportance(double al)
	{
		double transparency = 1.0 - al / 255.0;
		return (1.0 - transparency * transparency);
	}

	void SetUpArrays() {
		for (int i = specials; i<netsize; i++) {
			double temp;
			temp = pow(i / 255.0, 1.0 / gamma_correction) * 255.0;
			temp = round(temp);
			biasvalues[i] = temp;

			network[i].b = network[i].g = network[i].r = biasvalue(i * 256 / netsize);

			/*  Sets alpha values at 0 for dark pixels. */
			if (i < 16)
				network[i].al = (i * 16);
			else
				network[i].al = BYTE_MAX;

			freq[i] = 1.0 / netsize;
		}
	}

	UINT unbiasvalue(double temp)
	{
		if (temp < 0)
			return 0;

		temp = pow(temp / 255.0, gamma_correction) * BYTE_MAX;
		temp = floor((temp / 255.0 * 256.0));

		if (temp > BYTE_MAX)
			return BYTE_MAX;
		return (UINT)temp;
	}

	inline UINT round_biased(double temp)
	{
		if (temp < 0)
			return 0;
		temp = floor((temp / 255.0 * 256.0));

		if (temp > BYTE_MAX)
			return BYTE_MAX;
		return (UINT)temp;
	}

	void Altersingle(double alpha, UINT i, byte al, byte r, byte g, byte b) {
		double colorimp = 1.0;//0.5;// + 0.7 * colorimportance(al);

		alpha /= initalpha;

		/* alter hit neuron */
		network[i].al -= alpha * (network[i].al - al);
		network[i].b -= colorimp * alpha * (network[i].b - b);
		network[i].g -= colorimp * alpha * (network[i].g - g);
		network[i].r -= colorimp * alpha * (network[i].r - r);
	}

	void Alterneigh(UINT rad, UINT i, byte al, byte r, byte g, byte b) {
		int lo = i - rad;
		if (lo < 0)
			lo = 0;
		UINT hi = i + rad;
		if (hi > netsize - 1)
			hi = netsize - 1;

		UINT j = i + 1;
		int k = i - 1;
		double* q = radpower;
		while ((j <= hi) || (k >= lo)) {
			double a = (*(++q)) / alpharadbias;
			if (j <= hi) {
				network[j].al -= a * (network[j].al - al);
				network[j].b -= a * (network[j].b - b);
				network[j].g -= a * (network[j].g - g);
				network[j].r -= a * (network[j].r - r);
				j++;
			}
			if (k >= lo) {
				network[k].al -= a * (network[k].al - al);
				network[k].b -= a * (network[k].b - b);
				network[k].g -= a * (network[k].g - g);
				network[k].r -= a * (network[k].r - r);
				k--;
			}
		}
	}

	int Contest(byte al, byte r, byte g, byte b) {
		/* finds closest neuron (min dist) and updates freq */
		/* finds best neuron (min dist-bias) and returns position */
		/* for frequently chosen neurons, freq[i] is high and bias[i] is negative */
		/* bias[i] = gamma*((1/netsize)-freq[i]) */

		UINT i;
		double dist, a, betafreq;
		UINT bestpos, bestbiaspos;
		double bestd, bestbiasd;

		bestd = 1 << 30;
		bestbiasd = bestd;
		bestpos = 0;
		bestbiaspos = bestpos;

		/* Using colorimportance(al) here was causing problems with images that were close to monocolor.
		See bug reports: 3149791, 2938728, 2896731 and 2938710
		*/
		double colimp = 1.0; //colorimportance(al); 

		for (i = 0; i < netsize; i++) {
			double bestbiasd_biased = bestbiasd + bias[i];

			a = network[i].b - b;
			dist = abs(a) * colimp;
			a = network[i].r - r;
			dist += abs(a) * colimp;

			if (dist < bestd || dist < bestbiasd_biased) {
				a = network[i].g - g;
				dist += abs(a) * colimp;
				a = network[i].al - al;
				dist += abs(a);

				if (dist<bestd) {
					bestd = dist;
					bestpos = i;
				}
				if (dist < bestbiasd_biased) {
					bestbiasd = dist - bias[i];
					bestbiaspos = i;
				}
			}
			betafreq = freq[i] / (1 << betashift);
			freq[i] -= betafreq;
			bias[i] += betafreq * (1 << gammashift);
		}
		freq[bestpos] += beta;
		bias[bestpos] -= betagamma;
		return bestbiaspos;
	}

	void Learn(int samplefac) {
		UINT i, j;
		UINT rad, delta, samplepixels, stepIndex = 0;
		double radius, alpha;

		int pos = 0;
		int alphadec = 30 + ((samplefac - 1) / 3);
		int lengthcount = pixels.size();
		samplepixels = lengthcount / samplefac;
		delta = samplepixels / ncycles;  /* here's a problem with small images: samplepixels < ncycles => delta = 0 */
		if (delta == 0)
			delta = 1;        /* kludge to fix */
		alpha = initalpha;
		radius = initradius;

		rad = (UINT)radius;
		if (rad <= 1) rad = 0;
		for (i = 0; i < rad; i++)
			radpower[i] = floor(alpha * (((rad*rad - i * i) * radiusbias) / (rad * rad)));

		i = 0;
		while (i < samplepixels) {
			Color c(pixels[pos]);

			byte al = c.GetA();
			if (al) {
				byte b = c.GetB();
				byte g = c.GetG();
				byte r = c.GetR();

				j = Contest(al, r, g, b);

				Altersingle(alpha, j, al, r, g, b);
				if (rad)
					Alterneigh(rad, j, al, r, g, b);   /* alter neighbours */
			}

			pos += primes[stepIndex++ % 4];
			while (pos >= lengthcount)
				pos -= lengthcount;

			i++;
			if (i % delta == 0) {                    /* FPE here if delta=0*/
				alpha -= alpha / (double)alphadec;
				radius -= radius / (double)radiusdec;
				rad = (UINT)radius;
				if (rad <= 1) rad = 0;
				for (j = 0; j<rad; j++)
					radpower[j] = floor(alpha * (((rad * rad - j * j) * radiusbias) / (rad * rad)));
			}
		}
	}

	void Inxbuild(ColorPalette* pPalette) {
		int k = 0;
		if(hasTransparency)
			++k;
		
		for (; k < netsize; k++) {
			auto alpha = round_biased(network[k].al);
			auto argb = Color::MakeARGB(alpha, biasvalue(unbiasvalue(network[k].r)), biasvalue(unbiasvalue(network[k].g)), biasvalue(unbiasvalue(network[k].b)));
			pPalette->Entries[k] = argb;
		}
		if(hasTransparency)
			pPalette->Entries[0] = m_transparentColor;

		int previouscol = 0;
		int startpos = 0;

		for (int i = 0; i<netsize; i++) {		
			Color c(pPalette->Entries[i]);
			int smallpos = i;
			int smallval = c.GetG();			// index on g
											// find smallest in i..netsize-1
			for (int j = i + 1; j<netsize; j++) {
				Color c2(pPalette->Entries[j]);
				if (c2.GetG() < smallval) {		// index on g				
					smallpos = j;
					smallval = c2.GetG();	// index on g
				}
			}
			// swap p (i) and q (smallpos) entries
			if (i != smallpos)
				swap(pPalette->Entries[smallpos], pPalette->Entries[i]);

			// smallval entry is now in position i
			if (smallval != previouscol) {
				netindex[previouscol] = (startpos + i) >> 1;
				for (int j = previouscol + 1; j < smallval; j++)
					netindex[j] = i;
				previouscol = smallval;
				startpos = i;
			}
		}
		netindex[previouscol] = (startpos + maxnetpos) >> 1;
		for (int j = previouscol + 1; j<netsize; j++)
			netindex[j] = maxnetpos; // really 256
	}

	UINT Inxsearch(const ColorPalette* pPalette, ARGB argb) {
		UINT k = 0;
		Color c(argb);
		vector<short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			UINT nMaxColors = pPalette->Count;
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

	bool quantize_image(const vector<ARGB>& pixels, const ColorPalette* pPalette, UINT* qPixels, int width, int height, bool dither)
	{		
		if (dither) {
			bool odd_scanline = false;
			short *thisrowerr, *nextrowerr;
			int j, a_pix, r_pix, g_pix, b_pix, dir, two_val;
			const byte DJ = 4;
			const byte DITHER_MAX = 20;
			const int err_len = (width + 2) * DJ;
			byte range_tbl[DJ * 256] = { 0 };
			byte* range = &range_tbl[256];
			unique_ptr<short[]> erowErr = make_unique<short[]>(err_len);
			unique_ptr<short[]> orowErr = make_unique<short[]>(err_len);
			char dith_max_tbl[512] = { 0 };
			char* dith_max = &dith_max_tbl[256];
			short* erowerr = erowErr.get();
			short* orowerr = orowErr.get();

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
					qPixels[pixelIndex] = Inxsearch(pPalette, argb);

					Color c2(pPalette->Entries[qPixels[pixelIndex]]);
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
		}
		else {
			for (int i = 0; i < (width * height); i++)
				qPixels[i] = Inxsearch(pPalette, pixels[i]);
		}
		return true;
	}

	void NeuQuantizer::Clear() {
		memset((void*) network, 0, sizeof(network));
		hasTransparency = false;

		for (int i = 0; i < 32; ++i)
			radpower[i] = 0.0;
		for (int i = 0; i < netsize; ++i) {
			netindex[i] = 0;
			biasvalues[i] = 0;

			bias[i] = 0.0;
			freq[i] = 0.0;
		}

		pixels.clear();
		closestMap.clear();
	}

	bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, const UINT* qPixels)
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

	// The work horse for NeuralNet color quantizing.
	bool NeuQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		if (nMaxColors != 256)
			nMaxColors = 256;
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		bool r = true;
		int pixelIndex = 0;
		pixels.resize(bitmapWidth * bitmapHeight);
		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++) {
					Color color;
					pSource->GetPixel(x, y, &color);
					if(color.GetA() < BYTE_MAX)
						hasTransparency = true;
					if(color.GetA() == 0)
						m_transparentColor = color.GetValue();
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
					if (pixelAlpha < BYTE_MAX)
						hasTransparency = true;
					
					auto argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
					if (pixelAlpha == 0)
						m_transparentColor = argb;
					pixels[pixelIndex++] = argb;
				}

				pRowSource += strideSource;
			}

			pSource->UnlockBits(&data);
		}

		SetUpArrays();
		Learn(dither ? 5 : 1);

		auto pPaletteBytes = make_unique<byte[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		Inxbuild(pPalette);

		auto qPixels = make_unique<UINT[]>(bitmapWidth * bitmapHeight);
		quantize_image(pixels, pPalette, qPixels.get(), bitmapWidth, bitmapHeight, dither);
		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}
