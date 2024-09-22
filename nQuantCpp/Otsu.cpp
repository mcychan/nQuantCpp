/* Otsu's Image Segmentation Method
  Copyright (C) 2009 Tolga Birdal
  Copyright (c) 2018 - 2024 Miller Cy Chan
*/

#include "stdafx.h"
#include "Otsu.h"
#include "bitmapUtilities.h"
#include "GilbertCurve.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <unordered_map>

namespace OtsuThreshold
{
	BYTE alphaThreshold = 0xF;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, unsigned short > nearestMap;

	// function is used to compute the q values in the equation
	static float px(int init, int end, int* hist)
	{
		int sum = 0;
		for (int i = init; i <= end; ++i)
			sum += hist[i];

		return (float) sum;
	}

	// function is used to compute the mean values in the equation (mu)
	static float mx(int init, int end, int* hist)
	{
		int sum = 0;
		for (int i = init; i <= end; ++i)
			sum += i * hist[i];

		return (float) sum;
	}

	// finds the maximum element in a vector
	static short findMax(float* vec, int n)
	{
		float maxVec = 0;
		short idx = 0;

		for (int i = 1; i < n - 1; ++i) {
			if (vec[i] > maxVec) {
				maxVec = vec[i];
				idx = i;
			}
		}
		return idx;
	}

	// simply computes the image histogram
	void getHistogram(const vector<ARGB>& pixels, int* hist)
	{
		for (auto pixel : pixels) {
			Color c(pixel);
			if (c.GetA() <= alphaThreshold)
				continue;

			hist[c.GetR()]++;
			hist[c.GetG()]++;
			hist[c.GetB()]++;
		}
	}
		
	short getOtsuThreshold(const vector<ARGB>& pixels)
	{
		float vet[256] = { 0 };
		int hist[256] = { 0 };
		
		getHistogram(pixels, hist);

		// loop through all possible t values and maximize between class variance
		for (int k = 1; k != BYTE_MAX; ++k) {
			float p1 = px(0, k, hist);
			float p2 = px(k + 1, BYTE_MAX, hist);
			float p12 = p1 * p2;
			if (p12 == 0) 
				p12 = 1;
			float diff = (mx(0, k, hist) * p2) - (mx(k + 1, BYTE_MAX, hist) * p1);
			vet[k] = diff * diff / p12;
		}

		return findMax(vet, 256);
	}
	
	void threshold(const vector<ARGB>& pixels, vector<ARGB>& dest, short thresh, float weight = 1.0f)
	{
		auto maxThresh = (BYTE)thresh;
		if (thresh >= 200)
		{
			weight = .78f;
			maxThresh = (BYTE)(thresh * weight);
			thresh = 200;
		}

		auto minThresh = (BYTE)(thresh * weight);
		for (int i = 0; i < pixels.size(); ++i) {
			Color c(pixels[i]);
			if (m_transparentPixelIndex >= 0 && c.GetR() + c.GetG() + c.GetB() > maxThresh * 3)
				dest[i] = Color::MakeARGB(c.GetA(), BYTE_MAX, BYTE_MAX, BYTE_MAX);
			else if (m_transparentPixelIndex >= 0 || c.GetR() + c.GetG() + c.GetB() < minThresh * 3)
				dest[i] = Color::MakeARGB(c.GetA(), 0, 0, 0);
		}
	}

	vector<ARGB> cannyFilter(const UINT width, const vector<ARGB>& pixelsGray, double lowerThreshold, double higherThreshold) {
		const auto height = pixelsGray.size() / width;
		const auto area = (size_t)(width * height);

		vector<ARGB> pixelsCanny(area, Color::White);

		int gx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
		int gy[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
		auto G = make_unique<double[]>(area);
		vector<int> theta(area);
		auto largestG = 0.0;

		// perform canny edge detection on everything but the edges
		for (int i = 1; i < height - 1; ++i) {
			for (int j = 1; j < width - 1; ++j) {
				// find gx and gy for each pixel
				auto gxValue = 0.0;
				auto gyValue = 0.0;
				for (int x = -1; x <= 1; ++x) {
					for (int y = -1; y <= 1; ++y) {
						const Color c = pixelsGray[(i + x) * width +  j + y];
						gxValue += gx[1 - x][1 - y] * c.GetG();
						gyValue += gy[1 - x][1 - y] * c.GetG();
					}
				}

				const int center = i * width + j;
				// calculate G and theta
				G[center] = sqrt(pow(gxValue, 2) + pow(gyValue, 2));
				auto atanResult = atan2(gyValue, gxValue) * 180.0 / M_PI;
				theta[center] = (int)(180.0 + atanResult);

				if (G[center] > largestG)
					largestG = G[center];

				// setting the edges
				if (i == 1) {
					G[center - 1] = G[center];
					theta[center - 1] = theta[center];
				}
				else if (j == 1) {
					G[center - width] = G[center];
					theta[center - width] = theta[center];
				}
				else if (i == height - 1) {
					G[center + 1] = G[center];
					theta[center + 1] = theta[center];
				}
				else if (j == width - 1) {
					G[center + width] = G[center];
					theta[center + width] = theta[center];
				}

				// setting the corners
				if (i == 1 && j == 1) {
					G[center - width - 1] = G[center];
					theta[center - width - 1] = theta[center];
				}
				else if (i == 1 && j == width - 1) {
					G[center - width + 1] = G[center];
					theta[center - width + 1] = theta[center];
				}
				else if (i == height - 1 && j == 1) {
					G[center + width - 1] = G[center];
					theta[center + width - 1] = theta[center];
				}
				else if (i == height - 1 && j == width - 1) {
					G[center + width + 1] = G[center];
					theta[center + width + 1] = theta[center];
				}

				// to the nearest 45 degrees
				theta[center] = rint(theta[center] / 45) * 45;
			}
		}

		largestG *= .5;

		// non-maximum suppression
		for (int i = 1; i < height - 1; ++i) {
			for (int j = 1; j < width - 1; ++j) {
				const int center = i * width + j;
				if (theta[center] == 0 || theta[center] == 180) {
					if (G[center] < G[center - 1] || G[center] < G[center + 1])
						G[center] = 0;
				}
				else if (theta[center] == 45 || theta[center] == 225) {
					if (G[center] < G[center + width + 1] || G[center] < G[center - width - 1])
						G[center] = 0;
				}
				else if (theta[center] == 90 || theta[center] == 270) {
					if (G[center] < G[center + width] || G[center] < G[center - width])
						G[center] = 0;
				}
				else {
					if (G[center] < G[center + width - 1] || G[center] < G[center - width + 1])
						G[center] = 0;
				}

				auto grey = ~(BYTE)(G[center] * (255.0 / largestG));
				Color c(pixelsGray[center]);
				pixelsCanny[center] = Color::MakeARGB(c.GetA(), grey, grey, grey);
			}
		}

		int k = 0;
		auto minThreshold = lowerThreshold * largestG, maxThreshold = higherThreshold * largestG;
		do {
			for (int i = 1; i < height - 1; ++i) {
				for (int j = 1; j < width - 1; ++j) {
					const int center = i * width + j;
					if (G[center] < minThreshold)
						G[center] = 0;
					else if (G[center] >= maxThreshold)
						continue;
					else if (G[center] < maxThreshold) {
						G[center] = 0;
						for (int x = -1; x <= 1; ++x) {
							for (int y = -1; y <= 1; y++) {
								if (x == 0 && y == 0)
									continue;
								if (G[center + x * width + y] >= maxThreshold) {
									G[center] = higherThreshold * largestG;
									k = 0;
									x = 2;
									break;
								}
							}
						}
					}
					
					auto grey = ~(BYTE)(G[center] * 255.0 / largestG);
					Color c(pixelsGray[center]);
					pixelsCanny[center] = Color::MakeARGB(c.GetA(), grey, grey, grey);
				}
			}
		} while (k++ < 100);
		return pixelsCanny;
	}

	unsigned short nearestColorIndex(const ARGB* pPalette, const UINT nMaxColors, const ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return k;

		double mindist = INT_MAX;
		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette[i]);
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
		nearestMap[argb] = k;
		return k;
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	void convertToGrayScale(const vector<ARGB>& pixels, vector<ARGB>& dest)
	{
		float min1 = BYTE_MAX;
		float max1 = .0f;

		for (const auto& pixel : pixels)
		{
			int alfa = (pixel >> 24) & 0xff;
			if (alfa <= alphaThreshold)
				continue;

			int green = (pixel >> 8) & 0xff;
			if (min1 > green)
				min1 = green;

			if (max1 < green)
				max1 = green;
		}

		for (int i = 0; i < pixels.size(); ++i)
		{
			int alfa = (pixels[i] >> 24) & 0xff;
			if (alfa <= alphaThreshold)
				continue;

			int green = (pixels[i] >> 8) & 0xff;
			auto grey = (int)((green - min1) * (BYTE_MAX / (max1 - min1)));
			dest[i] = Color::MakeARGB(alfa, grey, grey, grey);
		}
	}

	Bitmap* Otsu::ConvertToGrayScale(Bitmap* pSrcImg)
	{
		auto iWidth = pSrcImg->GetWidth();
		auto iHeight = pSrcImg->GetHeight();

		auto pixelFormat = pSrcImg->GetPixelFormat();
		auto bitDepth = GetPixelFormatSize(pixelFormat);
		if (bitDepth != 32 && bitDepth != 24)
			pixelFormat = PixelFormat32bppARGB;

		Rect rect(0, 0, iWidth, iHeight);
		auto pSourceImg = pSrcImg->Clone(rect, pixelFormat);
		BitmapData data;
		Status status = pSourceImg->LockBits(&rect, ImageLockModeWrite, pSourceImg->GetPixelFormat(), &data);
		if (status != Ok)
			return pSrcImg;

		bitDepth = GetPixelFormatSize(pSourceImg->GetPixelFormat());
		auto DJ = (BYTE)(bitDepth >> 3);

		auto ptr = (LPBYTE)data.Scan0;

		float min1 = BYTE_MAX;
		float max1 = .0f;
		int remain = data.Stride - iWidth * DJ;

		for (int i = 0; i < iHeight; ++i)
		{
			for (int j = 0; j < iWidth; ++j)
			{
				if (DJ > 3 && ptr[3] <= alphaThreshold) {
					ptr += DJ;
					continue;
				}

				if (min1 > ptr[1])
					min1 = ptr[1];

				if (max1 < ptr[1])
					max1 = ptr[1];
				ptr += DJ;
			}
			ptr += remain;
		}

		ptr = (LPBYTE)data.Scan0;

		for (int i = 0; i < iHeight; ++i)
		{
			for (int j = 0; j < iWidth; ++j)
			{
				ptr[0] = ptr[1] = ptr[2] = (BYTE)((ptr[1] - min1) * (BYTE_MAX / (max1 - min1)));
				ptr += DJ;
			}
			ptr += remain;
		}

		pSourceImg->UnlockBits(&data);
		return pSourceImg;
	}


	bool Otsu::ConvertGrayScaleToBinary(Bitmap* pSrcImg, Bitmap* pDest, bool isGrayscale)
	{		
		const auto bitmapWidth = pSrcImg->GetWidth();
		const auto bitmapHeight = pSrcImg->GetHeight();
		const auto area = (size_t) (bitmapWidth * bitmapHeight);

		vector<ARGB> pixels(area);
		GrabPixels(pSrcImg, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor, alphaThreshold);

		auto pixelsGray = pixels;
		if (!isGrayscale)
			convertToGrayScale(pixels, pixelsGray);

		auto otsuThreshold = getOtsuThreshold(pixelsGray);
		auto lowerThreshold = 0.03, higherThreshold = 0.1;
		pixels = cannyFilter(bitmapWidth, pixelsGray, lowerThreshold, higherThreshold);
		threshold(pixelsGray, pixels, otsuThreshold);

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + 2 * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = 2;
		if (m_transparentPixelIndex >= 0) {
			pPalette->Entries[0] = m_transparentColor;
			pPalette->Entries[1] = Color::Black;
		}
		else {
			pPalette->Entries[0] = Color::Black;
			pPalette->Entries[1] = Color::White;
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette->Entries, pPalette->Count, nearestColorIndex, GetColorIndex, qPixels.get(), nullptr, 3.0f);
		if (m_transparentPixelIndex >= 0)
		{
			auto k = qPixels[m_transparentPixelIndex];
			if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}

		nearestMap.clear();
		pDest->SetPalette(pPalette);
		return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex >= 0);
	}
}
