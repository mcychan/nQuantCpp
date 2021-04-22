#pragma once

#include "stdafx.h"
#include "NeuQuantizer.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include <unordered_map>

namespace NeuralNet
{
	//====================
	// NeuralNet Color quantizing

	/* NeuQuant Neural-Net Quantization Algorithm
	* ------------------------------------------
	*
	* Copyright (c) 1994 Anthony Dekker
	* Copyright (c) 2019 Miller Cy Chan
	*
	* NEUQUANT Neural-Net quantization algorithm by Anthony Dekker, 1994.
	* See "Kohonen neural networks for optimal colour quantization"
	* in "Network: Computation in Neural Systems" Vol. 5 (1994) pp 351-367.
	* for research paper of the algorithm.
	* See also  https://www.researchgate.net/publication/232079905_Kohonen_neural_networks_for_optimal_colour_quantization
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

	double PR = .2126, PG = .7152, PB = .0722;
	const short specials = 3;		// number of reserved colours used
	const int ncycles = 115;			// no. of learning cycles
	const int radiusbiasshift = 8;
	const int radiusbias = 1 << radiusbiasshift;

	int netsize = 256;		// number of colours used	
	int maxnetpos = netsize - 1;
	int initrad = netsize >> 3;   // for 256 cols, radius starts at 32
	double initradius = initrad * 1.0;
	const int radiusdec = 30; // factor of 1/30 each cycle

	const short normal_learning_extension_factor = 2; /* normally learn twice as long */
	const short extra_long_colour_threshold = 40; /* learn even longer when under 40 */
	const short extra_long_divisor = 8; /* fraction of netsize for extra learning */

	const int alphabiasshift = 10;			// alpha starts at 1
	const int initalpha = 1 << alphabiasshift; // biased by 10 bits
	const int alpharadbshift = alphabiasshift + radiusbiasshift;
	const double alpharadbias = (double)(1 << alpharadbshift);
	const double exclusion_threshold = 0.5;

	const short REPEL_THRESHOLD = 16;          /* See repel_coincident()... */
	const short REPEL_STEP_DOWN = 1;              /* ... for an explanation of... */
	const short REPEL_STEP_UP = 4;                 /* ... how these points work. */
	unique_ptr<unsigned short[]> repel_points;

	/* defs for freq and bias */
	const int gammashift = 10;                  /* gamma = 1024 */
	const double gamma = (double)(1 << gammashift);
	const int betashift = 10;
	const double beta = (1.0 / (double)(1 << betashift));/* beta = 1/1024 */
	const double betagamma = (double)(1 << (gammashift - betashift));

	struct nq_pixel
	{
		double al, L, A, B;
	};

	unique_ptr<nq_pixel[]> network; // the network itself

	unique_ptr<unsigned short[]> netindex; // for network lookup - really 256

	unique_ptr<double[]> bias;  // bias and freq arrays for learning
	unique_ptr<double[]> freq;
	unique_ptr<double[]> radpower;

	double gamma_correction = 1.0;         // 1.0/2.2 usually

	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
	unordered_map<ARGB, unsigned short> nearestMap;

	inline double colorimportance(double al)
	{
		double transparency = 1.0 - al / 255.0;
		return (1.0 - transparency * transparency);
	}

	/* posneg(x) returns +1 if x is positive or zero, or -1 otherwise.
	*
	* XXX Posneg could be turned into a vector function in various ways.
	*/
	inline double posneg(double x) {
		if (x < 0)
			return -1.0;		
		return 1.0;
	}

	void SetUpArrays() {
		network = make_unique<nq_pixel[]>(netsize);
		netindex = make_unique<unsigned short[]>(max(netsize, 256));
		repel_points = make_unique<unsigned short[]>(max(netsize, 256));
		bias = make_unique<double[]>(netsize);
		freq = make_unique<double[]>(netsize);
		radpower = make_unique<double[]>(initrad);

		for (int i = specials; i < netsize; ++i) {
			network[i].L = network[i].A = network[i].B = i / netsize;

			/*  Sets alpha values at 0 for dark pixels. */
			if (i < 16)
				network[i].al = i * 16.0;
			else
				network[i].al = BYTE_MAX;

			freq[i] = 1.0 / netsize;
		}
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

	inline UINT round_biased(double temp)
	{
		if (temp < 0)
			return 0;
		temp = floor((temp / 255.0 * 256.0));

		if (temp > BYTE_MAX)
			return BYTE_MAX;
		return (UINT)temp;
	}

	void Altersingle(double alpha, UINT i, BYTE al, double L, double A, double B) {
		double colorimp = 1.0;//0.5;// + 0.7 * colorimportance(al);

		alpha /= initalpha;

		/* alter hit neuron */
		network[i].al -= alpha * (network[i].al - al);
		network[i].L -= colorimp * alpha * (network[i].L - L);
		network[i].A -= colorimp * alpha * (network[i].A - A);
		network[i].B -= colorimp * alpha * (network[i].B - B);
	}

	void Alterneigh(UINT rad, UINT i, BYTE al, double L, double A, double B) {
		int lo = i - rad;
		if (lo < 0)
			lo = 0;
		UINT hi = i + rad;
		if (hi > maxnetpos)
			hi = maxnetpos;

		UINT j = i + 1;
		int k = i - 1;
		auto learning_rate_table = radpower.get();
		while ((j <= hi) || (k >= lo)) {
			double learning_rate = (*(++learning_rate_table)) / alpharadbias;
			if (j <= hi) {
				network[j].al -= learning_rate * (network[j].al - al);
				network[j].L -= learning_rate * (network[j].L - L);
				network[j].A -= learning_rate * (network[j].A - A);
				network[j].B -= learning_rate * (network[j].B - B);
				j++;
			}
			if (k >= lo) {
				network[k].al -= learning_rate * (network[k].al - al);
				network[k].L -= learning_rate * (network[k].L - L);
				network[k].A -= learning_rate * (network[k].A - A);
				network[k].B -= learning_rate * (network[k].B - B);
				k--;
			}
		}
	}

	/* Repelcoincident(i) traverses the entire neural network, identifies any neurons that are close in colour to neuron i, and
	 * moves them away from neuron i.  Our definition of 'close' is being within exclusion_threshold of neuron i in each component.
	 *
	 * Because Repelcoincident() can be cpu-intensive, each neuron is given repel points (in the repel_points array) to limit how
	 * often a full repel pass is made.  Each time neuron i gets a full repel pass, it gets REPEL_STEP_UP points.  When the number
	 * of points it has exceeds REPEL_THRESHOLD, we stop doing full repel passes, and deduct REPEL_STEP_DOWN points instead.
	 * Eventually the number of repel points will eventually oscillate around the threshold.  With current settings, that means
	 * that only every 4th function call will result in a full pass.
 	*/
	void Repelcoincident(int i) {
		/* Use brute force to precompute the distance vectors between our neuron and each neuron. */

		if (repel_points[i] > REPEL_THRESHOLD) {
			repel_points[i] -= REPEL_STEP_DOWN;
			return;
		}

		auto diffs = make_unique<nq_pixel[]>(netsize);
		for (int vdx = 0; vdx < netsize; ++vdx) {
			diffs[vdx].al = abs(network[vdx].al - network[i].al);
			diffs[vdx].L = abs(network[vdx].L - network[i].L);
			diffs[vdx].A = abs(network[vdx].A - network[i].A);
			diffs[vdx].B = abs(network[vdx].B - network[i].B);
		}


		/* Identify which neurons are too close to neuron[i] and shift them away.
		 * */

		 /* repel_step is the amount we move the neurons away by in each component.
		  * radpower[0]/alpharadbias is similar to the proportion alterneigh uses.
		  * */
		double repel_step = exclusion_threshold * (radpower[0] / alpharadbias);

		for (int j = 0; j < netsize; ++j) {

			if (diffs[j].al < exclusion_threshold
				&& diffs[j].L < exclusion_threshold
				&& diffs[j].A < exclusion_threshold
				&& diffs[j].B < exclusion_threshold
				) {

				nq_pixel repel_vec;
				repel_vec.al = posneg(diffs[j].al), repel_vec.L = posneg(diffs[j].L), repel_vec.A = posneg(diffs[j].A), repel_vec.B = posneg(diffs[j].B);

				network[j].al += repel_vec.al * repel_step;
				network[j].L += repel_vec.L * repel_step;
				network[j].A += repel_vec.A * repel_step;
				network[j].B += repel_vec.B * repel_step;
			}
		}

		repel_points[i] += REPEL_STEP_UP;
	}

	int Contest(BYTE al, double L, double A, double B) {
		/* Calculate the component-wise differences between target_pix colour and every colour in the network, and weight according
		* to component relevance.
		*/
		auto diffs = make_unique<nq_pixel[]>(netsize);
		for (int vdx = 0; vdx < netsize; ++vdx) {
			diffs[vdx].al = abs(network[vdx].al - al);
			diffs[vdx].L = abs(network[vdx].L - L);
			diffs[vdx].A = abs(network[vdx].A - A);
			diffs[vdx].B = abs(network[vdx].B - B);
		}

		/* finds closest neuron (min dist) and updates freq */
		/* finds best neuron (min dist-bias) and returns position */
		/* for frequently chosen neurons, freq[i] is high and bias[i] is negative */
		/* bias[i] = gamma*((1/netsize)-freq[i]) */

		int bestpos = 0, bestbiaspos = bestpos;
		double bestd = INT_MAX, bestbiasd = bestd;

		/* Using colorimportance(al) here was causing problems with images that were close to monocolor.
		See bug reports: 3149791, 2938728, 2896731 and 2938710
		*/
		bool perfect = false; /* Is bestpos a perfect match for the colour. */

		for (int i = 0; i < netsize; ++i) {
			if (!perfect) {
				double bestbiasd_biased = bestbiasd + bias[i];
				double dist = 0;

				/*Detect perfect matches. */
				if (diffs[i].al < exclusion_threshold
					&& diffs[i].L < exclusion_threshold
					&& diffs[i].A < exclusion_threshold
					&& diffs[i].B < exclusion_threshold
					) {
					perfect = true;

					/* If we don't have a perfect match, the distance between the
						* target and the neuron is the manhattan distance. */
				}
				else {
					dist = diffs[i].L;
					dist += diffs[i].A;

					if (dist < bestd || dist < bestbiasd_biased) {
						dist += diffs[i].B;
						dist += diffs[i].al;						
					}
				}

				/* See if the current neuron is better. */
				if (dist < bestd) {
					bestd = dist;
					bestpos = i;
				}
				if (dist < bestbiasd_biased) {
					bestbiasd = dist - bias[i];
					bestbiaspos = i;
				}
			}

			/* Age (decay) the current neurons bias and freq values. */
			double betafreq = freq[i] * beta;
			freq[i] -= betafreq;
			bias[i] += betafreq * gamma;
		}

		/* Increase the freq and bias values for the chosen neuron. */
		freq[bestpos] += beta;
		bias[bestpos] -= betagamma;
		
		/* If our bestpos pixel is a 'perfect' match, we return bestpos, not bestbiaspos.  That is, we only decide to look at
		* bestbiaspos if the current target pixel wasn't a good enough match with the bestpos neuron, and there is some hope that
		* we can train the bestbiaspos neuron to become a better match. */
		if (perfect)
			return -bestpos - 1;  /* flag this was a perfect match */

		return bestbiaspos;
	}

	void Learn(const int samplefac, const vector<ARGB>& pixels) {
		UINT stepIndex = 0;

		int pos = 0;
		int alphadec = 30 + ((samplefac - 1) / 3);
		const UINT lengthcount = pixels.size();
		UINT samplepixels = lengthcount / samplefac;
		UINT delta = samplepixels / ncycles;  /* here's a problem with small images: samplepixels < ncycles => delta = 0 */
		if (delta == 0)
			delta = 1;        /* kludge to fix */
		double alpha = initalpha;
		double radius = initradius;
		if (netsize <= 2 * initradius)
			radius = netsize / 2;

		UINT rad = (UINT) radius;
		if (rad <= 1)
			rad = 0;

		for (UINT i = 0; i < rad; ++i)
			radpower[i] = floor(alpha * (((sqr(rad) - sqr(i)) * radiusbias) / sqr(rad)));

		UINT step = ((float)rand() / (float)RAND_MAX) * lengthcount;

		int learning_extension = normal_learning_extension_factor;
		if (netsize < extra_long_colour_threshold)
			learning_extension = 2 + ((extra_long_colour_threshold - netsize) / extra_long_divisor);

		UINT i = 0;
		while (i < learning_extension * samplepixels) {
			Color c(pixels[pos]);

			BYTE al = c.GetA();
			CIELABConvertor::Lab lab1;
			getLab(c, lab1);

			auto j = Contest(al, lab1.L, lab1.A, lab1.B);

			/* Determine if the colour was a perfect match.  j contains a factor encoded boolean. Horrible code to extract it. */
			bool was_perfect = (j < 0);
			j = (j < 0 ? -(j + 1) : j);

			Altersingle(alpha, j, al, lab1.L, lab1.A, lab1.B);
			if (rad && !was_perfect)
				Alterneigh(rad, j, al, lab1.L, lab1.A, lab1.B);   /* alter neighbours */
			else if (rad && was_perfect)
				Repelcoincident(j);  /* repel neighbours in colour space */

			pos += step;
			while (pos >= lengthcount)
				pos -= lengthcount;

			if (++i % delta == 0) {                    /* FPE here if delta=0*/
				alpha -= alpha / (learning_extension * (double) alphadec);
				radius -= radius / (double)radiusdec;
				rad = (UINT)radius;
				if (rad <= 1)
					rad = 0;
				for (UINT j = 0; j < rad; ++j)
					radpower[j] = floor(alpha * (((sqr(rad) - sqr(j)) * radiusbias) / sqr(rad)));
			}
		}
	}

	void Inxbuild(ColorPalette* pPalette) {
		UINT nMaxColors = pPalette->Count;		

		int previouscol = 0;
		int startpos = 0;

		for (int i = 0; i < nMaxColors; ++i) {
			Color c(pPalette->Entries[i]);
			int smallpos = i;
			auto smallval = network[i].L;			// index on L
											// find smallest in i..netsize-1
			for (int j = i + 1; j < nMaxColors; ++j) {
				if (network[j].L < smallval) {		// index on L				
					smallpos = j;
					smallval = network[j].L;	// index on L
				}
			}
			// swap p (i) and q (smallpos) entries
			if (i != smallpos)
				swap(network[smallpos], network[i]);

			// smallval entry is now in position i
			if (smallval != previouscol) {
				netindex[previouscol] = (startpos + i) >> 1;
				for (int j = previouscol + 1; j < smallval; ++j)
					netindex[j] = i;
				previouscol = smallval;
				startpos = i;
			}
		}

		netindex[previouscol] = (startpos + maxnetpos) >> 1;
		for (int j = previouscol + 1; j < netsize; ++j)
			netindex[j] = maxnetpos;

		for (UINT k = 0; k < nMaxColors; ++k) {
			CIELABConvertor::Lab lab1;
			lab1.alpha = round_biased(network[k].al);
			lab1.L = network[k].L, lab1.A = network[k].A, lab1.B = network[k].B;
			pPalette->Entries[k] = CIELABConvertor::LAB2RGB(lab1);
		}
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;
		
		unsigned short k = 0;
		Color c(argb);

		double mindist = INT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette->Entries[i]);
			double curdist = sqr(c2.GetA() - c.GetA());
			if (curdist > mindist)
				continue;

			getLab(c2, lab2);
			if (nMaxColors > 32) {
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

					curdist += sqr(lab2.B - lab1.B) / 2.0;
				}
			}
			else {
				curdist += sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.B - lab1.B);				
			}
			
			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		nearestMap[argb] = k;
		return k;
	}

	bool quantize_image(const vector<ARGB>& pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither)
			return dither_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		UINT pixelIndex = 0;
		for (UINT j = 0; j < height; ++j) {
			for (UINT i = 0; i < width; ++i)
				qPixels[pixelIndex++] = nearestColorIndex(pPalette, nMaxColors, pixels[pixelIndex]);
		}

		return true;
	}

	void Clear() {
		network.reset();
		netindex.reset();
		bias.reset();
		freq.reset();
		radpower.reset();

		pixelMap.clear();
		nearestMap.clear();
	}

	// The work horse for NeuralNet color quantizing.
	bool NeuQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		netsize = nMaxColors;		// number of colours used
		maxnetpos = netsize - 1;
		initrad = netsize < 8 ? 1 : (netsize >> 3);
		initradius = initrad * 1.0;

		SetUpArrays();
		Learn(dither ? 5 : 1, pixels);
		Inxbuild(pPalette);

		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			dithering_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight);
			Clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		if (hasSemiTransparency || nMaxColors <= 32)
			PR = PG = PB = 1;

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		quantize_image(pixels, pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}

		Clear();
		return ProcessImagePixels(pDest, pPalette, qPixels.get(), m_transparentPixelIndex >= 0);
	}
}
