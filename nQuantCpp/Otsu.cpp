#pragma once
/* Otsu's Image Segmentation Method
  Copyright (C) 2009 Tolga Birdal
  Copyright (c) 2018 - 2021 Miller Cy Chan
*/

#include "stdafx.h"
#include "Otsu.h"
#include "bitmapUtilities.h"
#include "GilbertCurve.h"
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
		short idx= 0;

		for (short i = 1; i < n - 1; ++i) {
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
	
	void threshold(vector<ARGB>& pixels, short thresh, float weight = 1.0f)
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
			if (c.GetR() + c.GetG() + c.GetB() > maxThresh * 3)
				pixels[i] = Color::MakeARGB(c.GetA(), BYTE_MAX, BYTE_MAX, BYTE_MAX);
			else if (m_transparentPixelIndex >= 0 || c.GetR() + c.GetG() + c.GetB() < minThresh * 3)
				pixels[i] = Color::MakeARGB(c.GetA(), 0, 0, 0);
		}
	}

	unsigned short nearestColorIndex(const ColorPalette* pPalette, const ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return k;

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
		nearestMap[argb] = k;
		return k;
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	void convertToGrayScale(vector<ARGB>& pixels)
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
			pixels[i] = Color::MakeARGB(alfa, grey, grey, grey);
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
		auto bitmapWidth = pSrcImg->GetWidth();
		auto bitmapHeight = pSrcImg->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSrcImg, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor, alphaThreshold);

		if (!isGrayscale)
			convertToGrayScale(pixels);

		auto otsuThreshold = getOtsuThreshold(pixels);
		threshold(pixels, otsuThreshold);

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
		Peano::GilbertCurve::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nearestColorIndex, GetColorIndex, qPixels.get());
		if (m_transparentPixelIndex >= 0)
		{
			auto k = qPixels[m_transparentPixelIndex];
			if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}

		nearestMap.clear();
		return ProcessImagePixels(pDest, pPalette, qPixels.get(), m_transparentPixelIndex >= 0);
	}
}
