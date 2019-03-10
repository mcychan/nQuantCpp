#pragma once
/* Multiobjective Image Color Quantization Algorithm Based on Self-Adaptive Hybrid Differential Evolution
Copyright (C) 2014-2016 Zhongbo Hu, Qinghua Su, Xuewen Xia
Copyright (c) 2018 Miller Cy Chan
* Adaptive algorithm k-means idea-accelerated differential evolution algorithm: each individual is a set of cluster centers, according to the probability K_probability
* The mean center of the K_number sub-cluster of this class replaces the original individual (following the acceptance criteria of parent-child competition)
* Iterate the differential evolution algorithm
* Initialize the individual is the pixel in the picture
* The program with the smallest error within the class + Maximum distance between classes + Minimize MSE */

#include "stdafx.h"
#include "MoDEQuantizer.h"
#include <map>

namespace MoDEQuant
{
	const double a1 = 0.1, a2 = 0.1, a3 = 0.8;  // Linear combination parameters
	const unsigned short K_number = 10;    // Number of cluster iterations
	const double K_probability = 0.01; // Probability of cluster iteration for each individual
	const unsigned short N = 100;           //  population size
	const double Max_diff = 100.0; // Maximum variation
	const byte SIDE = 3;
	const byte high = BYTE_MAX;
	const byte low = 0;        // the initial bounding region
	const int my_gens = 200;   //the generation number
	const int LOOP = 10;        //loop number
	const int seed[] = { 20436,18352,10994,26845,24435,29789,28299,11375,10222,9885,25855,4282,22102,29385,16014,32018,3200,11252,6227,5939,8712,12504,25965,6101,30359,1295,29533,19841,14690,2695,3503,16802,18931,28464,1245,13279,5676,8951,7280,24488,6537,27128,9320,16399,24997,24303,16862,17882,15360,31216 };

	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	map<ARGB, vector<short> > closestMap;

	inline const int getARGBIndex(const Color& c)
	{
		if (hasSemiTransparency)
			return (c.GetA() & 0xF0) << 8 | (c.GetR() & 0xF0) << 4 | (c.GetG() & 0xF0) | (c.GetB() >> 4);
		if (m_transparentPixelIndex >= 0)
			return (c.GetA() & 0x80) << 8 | (c.GetR() & 0xF8) << 7 | (c.GetG() & 0xF8) << 2 | (c.GetB() >> 3);
		return (c.GetR() & 0xF8) << 8 | (c.GetG() & 0xFC) << 3 | (c.GetB() >> 3);
	}

	inline double sqr(double value)
	{
		return value * value;
	}

	inline double rand1()
	{
		return (double)rand() / (RAND_MAX + 1.0);
	}

	unsigned short find_nn(const vector<double>& data, const Color& c)
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		unsigned short temp_k = nMaxColors;  //Record the ith pixel is divided into classes in the center of the temp_k
		double idis = INT_MAX;
		int nIndex = 0;
		for (unsigned short k = 0; k<nMaxColors; ++k) {
			double iid0 = sqr(data[nIndex] - c.GetR());
			double iid1 = sqr(data[nIndex + 1] - c.GetG());
			double iid2 = sqr(data[nIndex + 2] - c.GetB());
			double iidis = iid0 + iid1 + iid2;
			if (iidis < idis) {
				idis = iidis;
				temp_k = k;   //Record the ith pixel is divided into classes in the center of the temp_k
			}
			nIndex += SIDE;
		}
		return temp_k;
	}

	void updateCentroids(vector<double>& data, double* temp_x, const int* temp_x_number)
	{
		const unsigned short nMaxColors = data.size() / SIDE;

		for (unsigned short i = 0; i<nMaxColors; ++i) { //update classes and centroids
			if (temp_x_number[i] > 0) {
				data[SIDE * i] = temp_x[SIDE * i] / temp_x_number[i];
				data[SIDE * i + 1] = temp_x[SIDE * i + 1] / temp_x_number[i];
				data[SIDE * i + 2] = temp_x[SIDE * i + 2] / temp_x_number[i];
			}
		}
	}

	double evaluate1(const vector<ARGB>& pixels, const vector<double>& data)
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		double dis_sum = 0.0;
		UINT nSize = pixels.size();
		auto k_class_Number = make_unique<int[]>(nMaxColors); //store distance of each class and related parameters
		auto k_class_dis = make_unique<double[]>(nMaxColors);
		for (UINT i = 0; i < nSize; ++i) {
			Color c(pixels[i]);
			auto k_class_Temp = find_nn(data, c);
			if (k_class_Temp < nMaxColors) {
				double id0 = sqr(data[SIDE * k_class_Temp] - c.GetR());
				double id1 = sqr(data[SIDE * k_class_Temp + 1] - c.GetG());
				double id2 = sqr(data[SIDE * k_class_Temp + 2] - c.GetB());

				k_class_dis[k_class_Temp] += (id0 + id1 + id2);
				k_class_Number[k_class_Temp]++;
			}
		}
		for (unsigned short j = 0; j < nMaxColors; ++j) {
			double TT = k_class_dis[j] / k_class_Number[j];
			if (dis_sum < TT)
				dis_sum = TT;
		}
		return dis_sum;
	}

	// Adaptation function designed for multiple targets (1): the minimum value of each inner class distance is the smallest
	double evaluate1_K(const vector<ARGB>& pixels, vector<double>& data)  //Adaptive value function with K-means variation
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		UINT nSize = pixels.size();
		auto temp_i_k = make_unique<int[]>(nSize);

		for (int ii = 0; ii<K_number; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(SIDE * nMaxColors);  //store average value centroid
			for (unsigned short i = 0; i<nSize; i++) {
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c);
				if (temp_k < nMaxColors) {
					temp_x_number[temp_k]++;
					temp_x[SIDE * temp_k] += c.GetR();     //Put each pixel of the original image into categories and put them in an array
					temp_x[SIDE * temp_k + 1] += c.GetG();
					temp_x[SIDE * temp_k + 2] += c.GetB();
					temp_i_k[i] = temp_k;
				}
			}

			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		auto k_class_dis = make_unique<double[]>(nMaxColors);;
		auto k_class_Number = make_unique<int[]>(nMaxColors);;
		for (UINT i = 0; i<nSize; ++i) {
			Color c(pixels[i]);
			int j = SIDE * temp_i_k[i];
			int jj = temp_i_k[i];
			k_class_Number[jj]++;
			double iid0 = sqr(data[j] - c.GetR());
			double iid1 = sqr(data[j + 1] - c.GetG());
			double iid2 = sqr(data[j + 2] - c.GetB());
			k_class_dis[jj] += (iid0 + iid1 + iid2);
		}

		double dis_sum = 0.0;
		for (unsigned short j = 0; j < nMaxColors; ++j) {
			double TT = k_class_dis[j] / k_class_Number[j];
			if (dis_sum < TT)
				dis_sum = TT;
		}
		return dis_sum;
	}

	double evaluate2(const vector<ARGB>& pixels, vector<double>& data, const int K_num = 1)  //Adaptive value function with K-means variation
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		UINT nSize = pixels.size();

		for (int ii = 0; ii<K_num; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(SIDE * nMaxColors);  //store average value centroid
			for (UINT i = 0; i<nSize; ++i) {
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c);
				if (temp_k < nMaxColors) {
					temp_x_number[temp_k]++;
					temp_x[SIDE * temp_k] += c.GetR();     //Put each pixel of the original image into categories and put them in an array
					temp_x[SIDE * temp_k + 1] += c.GetG();
					temp_x[SIDE * temp_k + 2] += c.GetB();
				}
			}

			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		double Temp_dis = INT_MAX;
		for (unsigned short i = 0; i<nMaxColors - 1; ++i) {     // Calculate the fitness value after several clusters（class with minimum distance）
			for (unsigned short j = i + 1; j<nMaxColors; ++j) {
				double iid0 = sqr(data[SIDE * i] - data[SIDE * j]);
				double iid1 = sqr(data[SIDE * i + 1] - data[SIDE * j + 1]);
				double iid2 = sqr(data[SIDE * i + 2] - data[SIDE * j + 2]);
				double T_Temp_dis = iid0 + iid1 + iid2;
				if (Temp_dis > T_Temp_dis)
					Temp_dis = T_Temp_dis;
			}
		}
		return (Temp_dis);
	}

	// designed for multi objective application function(2)：to maximize the minimum distance of class
	double evaluate2_K(const vector<ARGB>& pixels, vector<double>& data)  //Adaptive value function with K-means variation
	{
		return evaluate2(pixels, data, K_number);
	}

	//designed for multi objective application function(3) MSE
	double evaluate3(const vector<ARGB>& pixels, const vector<double>& data)
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		double dis_sum = 0.0;
		UINT nSize = pixels.size();
		for (UINT i = 0; i < nSize; ++i) {
			Color c(pixels[i]);
			auto k_class_Temp = find_nn(data, c);
			if (k_class_Temp < nMaxColors) {
				double id0 = sqr(data[k_class_Temp * SIDE] - c.GetR());
				double id1 = sqr(data[k_class_Temp * SIDE + 1] - c.GetG());
				double id2 = sqr(data[k_class_Temp * SIDE + 2] - c.GetB());
				dis_sum += (id0 + id1 + id2);
			}
		}

		return dis_sum / nSize;
	}

	double evaluate3_K(const vector<ARGB>& pixels, vector<double>& data)  //Adaptive value function with K-means variation
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		UINT nSize = pixels.size();
		auto temp_i_k = make_unique<int[]>(nSize);

		for (int ii = 0; ii<K_number; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(SIDE * nMaxColors);  //store average value centroid
			for (UINT i = 0; i<nSize; i++) {
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c);
				if (temp_k < nMaxColors) {
					temp_x_number[temp_k]++;
					temp_x[SIDE * temp_k] += c.GetR();     //Put each pixel of the original image into categories and put them in an array
					temp_x[SIDE * temp_k + 1] += c.GetG();
					temp_x[SIDE * temp_k + 2] += c.GetB();
					temp_i_k[i] = temp_k;
				}
			}

			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		double dis_sum = 0.0;
		for (UINT i = 0; i<nSize; ++i) {
			Color c(pixels[i]);
			int j = SIDE * temp_i_k[i];
			double iid0 = sqr(data[j] - c.GetR());
			double iid1 = sqr(data[j + 1] - c.GetG());
			double iid2 = sqr(data[j + 2] - c.GetB());
			dis_sum += (iid0 + iid1 + iid2);
		}

		return dis_sum / nSize;
	}

	int moDEquan(const vector<ARGB>& pixels, ColorPalette* pPalette, const unsigned short nMaxColors)
	{
		UINT nSizeInit = pixels.size();
		UINT D = nMaxColors * SIDE;
		vector<vector<double> > x1(N), x2(N);
		for (int i = 0; i<N; ++i) {
			x1[i].resize(D);
			x2[i].resize(D);
		}

		double cost[N];
		vector<double> bestx(D);

		for (int ii = 0; ii<LOOP; ++ii) { // loop n times to test
			double F = 0.5, CR = 0.6, BVATG = INT_MAX;
			srand(seed[ii]);
			for (int i = 0; i<N; ++i) {            //the initial population 
				for (UINT j = 0; j<D; ++j) {
					if (j % SIDE)
						continue;
					int TempInit = int(rand1() * nSizeInit);
					Color c(pixels[TempInit]);
					x1[i][j] = c.GetR();
					x1[i][j + 1] = c.GetG();
					x1[i][j + 2] = c.GetB();
				}

				cost[i] = a1 * evaluate1(pixels, x1[i]);
				cost[i] -= a2 * evaluate2(pixels, x1[i]);
				cost[i] += a3 * evaluate3(pixels, x1[i]) + 1000;

				if (cost[i] < BVATG) {
					BVATG = cost[i];
					for (UINT j = 0; j<D; ++j)
						bestx[j] = x1[i][j];
				}
			}

			for (int g = 1; g <= my_gens; ++g) { //generation loop
				for (int i = 0; i<N; ++i) {
					if ((rand1() < K_probability)) { // individual according to probability to perform clustering					

						double temp_costx1 = cost[i];
						cost[i] = a1 * evaluate1_K(pixels, x1[i]);
						cost[i] -= a2 * evaluate2_K(pixels, x1[i]);
						cost[i] += a3 * evaluate3_K(pixels, x1[i]) + 1000; // clustered and changed the original data of x1[i]

						if (cost[i] >= temp_costx1)
							cost[i] = temp_costx1;
						else if (BVATG >= cost[i]) {
							BVATG = cost[i];
							// update with the best individual
							for (UINT j = 0; j<D; ++j)
								bestx[j] = x1[i][j];
						}
					}
					else { // Differential Evolution
						int d, b;
						do {
							d = (int)(rand1() * N);
						} while (d == i);
						do {
							b = (int)(rand1() * N);
						} while (b == d || b == i);

						int jr = (int)(rand1() * D); // every individual update control parameters
						if (rand1() < 0.1) {
							F = 0.1 + rand1() * 0.9;
							CR = rand1();
						}

						for (UINT j = 0; j<D; ++j) {
							if (rand1() <= CR || j == jr) {
								double diff = (x1[d][j] - x1[b][j]);
								if (diff > Max_diff)
									diff -= Max_diff;
								if (diff > Max_diff)
									diff -= Max_diff;
								if (diff < -Max_diff)
									diff += Max_diff;
								if (diff < -Max_diff)
									diff += Max_diff;
								x2[i][j] = bestx[j] + F * diff;

								// periodic mode
								if (x2[i][j] < low)
									x2[i][j] = high - (low - x2[i][j]);
								else if (x2[i][j] > high)
									x2[i][j] = low + (x2[i][j] - high);
							}
							else
								x2[i][j] = x1[i][j];
						}

						double score = a1 * evaluate1(pixels, x1[i]);
						score -= a2 * evaluate2(pixels, x1[i]);
						if (score > cost[i])
							continue;
						score += a3 * evaluate3(pixels, x1[i]) + 1000;
						if (score > cost[i])
							continue;

						cost[i] = score;
						if (BVATG >= score) {
							BVATG = score;
							for (UINT j = 0; j < D; ++j)
								bestx[j] = x1[i][j] = x2[i][j];
						}
						else {
							for (UINT j = 0; j < D; ++j)
								x1[i][j] = x2[i][j];
						}
					}
				}

			}
		}

		/* Fill palette */
		UINT j = 0;
		for (UINT k = 0; k<nMaxColors; ++k) {
			pPalette->Entries[k] = Color::MakeARGB(BYTE_MAX, static_cast<byte>(bestx[j]), static_cast<byte>(bestx[j + 1]), static_cast<byte>(bestx[j + 2]));

			if (m_transparentPixelIndex >= 0 && Color(pPalette->Entries[k]).ToCOLORREF() == Color(m_transparentColor).ToCOLORREF()) {
				pPalette->Entries[k] = m_transparentColor;
				swap(pPalette->Entries[0], pPalette->Entries[k]);
			}
			j += SIDE;
		}

		return 0;
	}

	UINT nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		UINT k = 0;
		Color c(argb);

		UINT curdist, mindist = INT_MAX;
		for (short i = 0; i < nMaxColors; i++) {
			Color c2(pPalette->Entries[i]);
			int adist = sqr(c2.GetA() - c.GetA());
			curdist = adist;
			if (curdist > mindist)
				continue;

			int rdist = sqr(c2.GetR() - c.GetR());
			curdist += rdist;
			if (curdist > mindist)
				continue;

			int gdist = sqr(c2.GetG() - c.GetG());
			curdist += gdist;
			if (curdist > mindist)
				continue;

			int bdist = sqr(c2.GetB() - c.GetB());
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
			closest[2] = closest[3] = INT_MAX;

			for (; k < nMaxColors; k++) {
				Color c2(pPalette->Entries[k]);
				closest[4] = sqr(c.GetA() - c2.GetA()) + sqr(c.GetR() - c2.GetR()) + sqr(c.GetG() - c2.GetG()) + sqr(c.GetB() - c2.GetB());
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

			if (closest[3] == INT_MAX)
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
			cerr << "Cannot write image" << endl;
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
			cerr << "Cannot write image" << endl;
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

	bool MoDEQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
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

		if (nMaxColors > 2)
			moDEquan(pixels, pPalette, nMaxColors);
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette->Entries[0] = Color::Transparent;
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
			pPalette->Entries[k] = m_transparentColor;
		}
		closestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}