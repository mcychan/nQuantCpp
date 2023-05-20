/* Multiobjective Image Color Quantization Algorithm Based on Self-Adaptive Hybrid Differential Evolution
Copyright (C) 2014-2016 Zhongbo Hu, Qinghua Su, Xuewen Xia
Copyright (c) 2018 - 2021 Miller Cy Chan
* Adaptive algorithm k-means idea-accelerated differential evolution algorithm: each individual is a set of cluster centers, according to the probability K_probability
* The mean center of the K_number sub-cluster of this class replaces the original individual (following the acceptance criteria of parent-child competition)
* Iterate the differential evolution algorithm
* Initialize the individual is the pixel in the picture
* The program with the smallest error within the class + Maximum distance between classes + Minimize MSE */

#include "stdafx.h"
#include "MoDEQuantizer.h"
#include "bitmapUtilities.h"
#include <ctime>
#include <iomanip>      // std::setprecision
#include <unordered_map>

namespace MoDEQuant
{
	const double a1 = 0.1, a2 = 0.05, a3 = 0.01;  // Linear combination parameters
	const unsigned short K_number = 10;    // Number of cluster iterations
	const double K_probability = 0.05; // Probability of cluster iteration for each individual
	const unsigned short N = 100;           //  population size
	const double Max_diff = 100.0; // Maximum variation	
	const BYTE high = BYTE_MAX;
	const BYTE low = 0;        // the initial bounding region
	const int my_gens = 200;   //the generation number
	const int LOOP = 5;        //loop number
	const int seed[50] = { 20436,18352,10994,26845,24435,29789,28299,11375,10222,9885,25855,4282,22102,29385,16014,32018,3200,11252,6227,5939,8712,12504,25965,6101,30359,1295,29533,19841,14690,2695,3503,16802,18931,28464,1245,13279,5676,8951,7280,24488,6537,27128,9320,16399,24997,24303,16862,17882,15360,31216 };

	BYTE SIDE = 3;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, vector<unsigned short> > closestMap;

	inline double rand1()
	{
		return (double)rand() / (RAND_MAX + 1.0);
	}

	unsigned short find_nn(const vector<double>& data, const Color& c, double& idis)
	{
		auto argb = c.GetValue();
		const unsigned short nMaxColors = data.size() / SIDE;
		unsigned short temp_k = nMaxColors;  //Record the ith pixel is divided into classes in the center of the temp_k			
		for (unsigned short k = 0; k < nMaxColors; ++k) {
			double iidis = sqr(data[SIDE * k] - c.GetB());
			if (iidis >= idis)
				continue;

			iidis += sqr(data[SIDE * k + 1] - c.GetG());
			if (iidis >= idis)
				continue;

			iidis += sqr(data[SIDE * k + 2] - c.GetR());
			if (iidis >= idis)
				continue;

			if (hasSemiTransparency) {
				iidis += sqr(data[SIDE * k + 3] - c.GetA());
				if (iidis >= idis)
					continue;
			}

			idis = iidis;
			temp_k = k;   //Record the ith pixel is divided into classes in the center of the temp_k
		}
		return temp_k;
	}

	void updateCentroids(vector<double>& data, double* temp_x, const int* temp_x_number)
	{
		const unsigned short nMaxColors = data.size() / SIDE;

		for (unsigned short i = 0; i < nMaxColors; ++i) { //update classes and centroids
			if (temp_x_number[i] > 0) {
				data[SIDE * i] = temp_x[SIDE * i] / temp_x_number[i];
				data[SIDE * i + 1] = temp_x[SIDE * i + 1] / temp_x_number[i];
				data[SIDE * i + 2] = temp_x[SIDE * i + 2] / temp_x_number[i];
				if (hasSemiTransparency)
					data[SIDE * i + 3] = temp_x[SIDE * i + 3] / temp_x_number[i];
			}
		}
	}

	double evaluate1(const vector<ARGB>& pixels, const vector<double>& data)
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		UINT nSize = pixels.size();
		auto k_class_Number = make_unique<int[]>(nMaxColors); //store distance of each class and related parameters
		auto k_class_dis = make_unique<double[]>(nMaxColors);
		for (UINT i = 0; i < nSize; ++i) {
			double idis = INT_MAX;
			auto k_class_Temp = find_nn(data, pixels[i], idis);
			k_class_dis[k_class_Temp] += _sqrt(idis);
			k_class_Number[k_class_Temp]++;
		}

		double dis_sum = 0.0;
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

		for (int ii = 0; ii < K_number; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(data.size());  //store average value centroid
			for (UINT i = 0; i < nSize; ++i) {
				double idis = INT_MAX;
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c, idis);
				temp_x_number[temp_k]++;
				temp_x[SIDE * temp_k] += c.GetB();     //Put each pixel of the original image into categories and put them in an array
				temp_x[SIDE * temp_k + 1] += c.GetG();
				temp_x[SIDE * temp_k + 2] += c.GetR();
				if (hasSemiTransparency)
					temp_x[SIDE * temp_k + 3] += c.GetA();
				temp_i_k[i] = temp_k;
			}
			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		auto k_class_dis = make_unique<double[]>(nMaxColors);
		auto k_class_Number = make_unique<int[]>(nMaxColors);
		for (UINT i = 0; i < nSize; ++i) {
			Color c(pixels[i]);
			int j = SIDE * temp_i_k[i];
			int jj = temp_i_k[i];
			k_class_Number[jj]++;
			double iid0 = sqr(data[j] - c.GetB());
			double iid1 = sqr(data[j + 1] - c.GetG());
			double iid2 = sqr(data[j + 2] - c.GetR());
			double iid3 = hasSemiTransparency ? sqr(data[j + 3] - c.GetA()) : 0;
			k_class_dis[jj] += _sqrt(iid0 + iid1 + iid2 + iid3);
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

		for (int ii = 0; ii < K_num; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(data.size());  //store average value centroid
			for (UINT i = 0; i < nSize; ++i) {
				double idis = INT_MAX;
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c, idis);
				temp_x_number[temp_k]++;
				temp_x[SIDE * temp_k] += c.GetB();     //Put each pixel of the original image into categories and put them in an array
				temp_x[SIDE * temp_k + 1] += c.GetG();
				temp_x[SIDE * temp_k + 2] += c.GetR();
				if (hasSemiTransparency)
					temp_x[SIDE * temp_k + 3] += c.GetA();
			}
			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		double Temp_dis = INT_MAX;
		for (unsigned short i = 0; i < nMaxColors - 1; ++i) {     // Calculate the fitness value after several clusters（class with minimum distance）
			for (unsigned short j = i + 1; j < nMaxColors; ++j) {
				double T_Temp_dis = sqr(data[SIDE * i] - data[SIDE * j]);
				if (T_Temp_dis > Temp_dis)
					continue;

				T_Temp_dis += sqr(data[SIDE * i + 1] - data[SIDE * j + 1]);
				if (T_Temp_dis > Temp_dis)
					continue;

				T_Temp_dis += sqr(data[SIDE * i + 2] - data[SIDE * j + 2]);
				if (T_Temp_dis > Temp_dis)
					continue;

				if (hasSemiTransparency) {
					T_Temp_dis += sqr(data[SIDE * i + 3] - data[SIDE * j + 3]);
					if (T_Temp_dis > Temp_dis)
						continue;
				}

				Temp_dis = T_Temp_dis;
			}
		}
		return _sqrt(Temp_dis);
	}

	// designed for multi objective application function(2)：to maximize the minimum distance of class
	double evaluate2_K(const vector<ARGB>& pixels,  vector<double>& data)  //Adaptive value function with K-means variation
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
			double idis = INT_MAX;
			auto k_class_Temp = find_nn(data, pixels[i], idis);
			if (k_class_Temp < nMaxColors)
				dis_sum += _sqrt(idis);
		}

		return dis_sum / nSize;
	}

	double evaluate3_K(const vector<ARGB>& pixels, vector<double>& data)  //Adaptive value function with K-means variation
	{
		const unsigned short nMaxColors = data.size() / SIDE;
		UINT nSize = pixels.size();
		auto temp_i_k = make_unique<int[]>(nSize);

		for (int ii = 0; ii < K_number; ++ii) {
			auto temp_x_number = make_unique<int[]>(nMaxColors);  //store the pixel count of each class
			auto temp_x = make_unique<double[]>(data.size());  //store average value centroid
			for (UINT i = 0; i < nSize; ++i) {
				double idis = INT_MAX;
				Color c(pixels[i]);
				auto temp_k = find_nn(data, c, idis);
				temp_x_number[temp_k]++;
				temp_x[SIDE * temp_k] += c.GetB();     //Put each pixel of the original image into categories and put them in an array
				temp_x[SIDE * temp_k + 1] += c.GetG();
				temp_x[SIDE * temp_k + 2] += c.GetR();
				if (hasSemiTransparency)
					temp_x[SIDE * temp_k + 3] += c.GetA();
				temp_i_k[i] = temp_k;
			}
			updateCentroids(data, temp_x.get(), temp_x_number.get());
		}

		double dis_sum = 0.0;
		for (UINT i = 0; i < nSize; ++i) {
			Color c(pixels[i]);
			int j = SIDE * temp_i_k[i];
			double iid0 = sqr(data[j] - c.GetB());
			double iid1 = sqr(data[j + 1] - c.GetG());
			double iid2 = sqr(data[j + 2] - c.GetR());
			double iid3 = hasSemiTransparency ? sqr(data[j + 3] - c.GetA()) : 0;
			dis_sum += _sqrt(iid0 + iid1 + iid2 + iid3);
		}

		return dis_sum / nSize;
	}

	int moDEquan(const vector<ARGB>& pixels, ColorPalette* pPalette, const unsigned short nMaxColors)
	{
		const BYTE INCR_STEP = 1;
		const float INCR_PERC = INCR_STEP * 100.0f / my_gens;
		const clock_t begin = clock();
		cout << std::setprecision(1) << std::fixed;

		const UINT nSizeInit = pixels.size();
		const UINT D = nMaxColors * SIDE;
		auto x1 = make_unique<vector<double>[]>(N);
		auto x2 = make_unique<vector<double>[]>(N);
		for (int i = 0; i < N; ++i) {
			x1[i].resize(D);
			x2[i].resize(D);
		}

		double cost[N];
		auto bestx = make_unique<double[]>(D);

		float percCompleted = 0;
		const int ii = LOOP - 1;

		double F = 0.5, CR = 0.6, BVATG = INT_MAX;
		srand(seed[ii]);
		for (int i = 0; i < N; ++i) {            //the initial population 
			for (UINT j = 0; j < D; j += SIDE) {
				int TempInit = static_cast<int>(rand1() * nSizeInit);
				Color c(pixels[TempInit]);
				x1[i][j] = c.GetR();
			}

			cost[i] = a1 * evaluate1(pixels, x1[i]);
			cost[i] -= a2 * evaluate2(pixels, x1[i]);
			cost[i] += a3 * evaluate3(pixels, x1[i]) + 1000;

			if (cost[i] < BVATG) {
				BVATG = cost[i];
				for (UINT j = 0; j < D; ++j)
					bestx[j] = x1[i][j];
			}
		}

		for (int g = 0; g < my_gens; ++g) { //generation loop				
			if (g % INCR_STEP == 0) {
				int elapsed_secs = int(clock() - begin) / CLOCKS_PER_SEC;
				cout << "\rMultiobjective CQ ALGO Based on Self-Adaptive Hybrid DE: " << percCompleted << "% COMPL (" << elapsed_secs << " sec)" << std::flush;
				percCompleted += INCR_PERC;
			}

			for (int i = 0; i < N; ++i) {
				if (rand1() < K_probability) { // individual according to probability to perform clustering
					double temp_costx1 = cost[i];
					cost[i] = a1 * evaluate1_K(pixels, x1[i]);
					cost[i] -= a2 * evaluate2_K(pixels, x1[i]);
					cost[i] += a3 * evaluate3_K(pixels, x1[i]) + 1000; // clustered and changed the original data of x1[i]

					if (cost[i] >= temp_costx1)
						cost[i] = temp_costx1;
					else if (BVATG >= cost[i]) {
						BVATG = cost[i];
						// update with the best individual
						for (UINT j = 0; j < D; ++j)
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

					for (UINT j = 0; j < D; ++j) {
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

		int elapsed_secs = int(clock() - begin) / CLOCKS_PER_SEC;
		cout << "\rMultiobjective CQ ALGO Based on Self-Adaptive Hybrid DE: Well done!! (" << elapsed_secs << " sec)" << endl;

		/* Fill palette */
		UINT j = 0;
		for (unsigned short k = 0; k < nMaxColors; ++k, j += SIDE)
			pPalette->Entries[k] = Color::MakeARGB(hasSemiTransparency ? static_cast<BYTE>(bestx[j + 3]) : BYTE_MAX, static_cast<BYTE>(bestx[j + 2]), static_cast<BYTE>(bestx[j + 1]), static_cast<BYTE>(bestx[j]));

		return 0;
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, ARGB argb, const UINT pos)
	{
		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= 0)
			c = m_transparentColor;

		UINT mindist = INT_MAX;
		const auto nMaxColors = pPalette->Count;
		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette->Entries[i]);
			UINT curdist = sqr(c2.GetA() - c.GetA());
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
		UINT k = 0;
		Color c(argb);
		const auto nMaxColors = pPalette->Count;

		vector<unsigned short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = INT_MAX;

			for (; k < nMaxColors; ++k) {
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

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither)
			return dither_image(pixels, pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (int j = 0; j < height; ++j) {
			for (int i = 0; i < width; ++i)
				qPixels[pixelIndex++] = ditherFn(pPalette, pixels[pixelIndex], i + j);
		}

		return true;
	}

	bool MoDEQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		hasSemiTransparency = false;
		m_transparentPixelIndex = -1;
		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor, 0xF, nMaxColors);

		SIDE = hasSemiTransparency ? 4 : 3;
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
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

		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			dithering_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight);
			closestMap.clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		closestMap.clear();

		return ProcessImagePixels(pDest, pPalette, qPixels.get(), m_transparentPixelIndex >= 0);
	}

}