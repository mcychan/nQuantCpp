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
	BYTE alphaThreshold = 0;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, unsigned short > nearestMap;

	// function is used to compute the q values in the equation
	static float Px(int init, int end, int* hist)
	{
		int sum = 0;
		int i;
		for (i = init; i <= end; ++i)
			sum += hist[i];

		return (float) sum;
	}

	// function is used to compute the mean values in the equation (mu)
	static float Mx(int init, int end, int* hist)
	{
		int sum = 0;
		int i;
		for (i = init; i <= end; ++i)
			sum += i * hist[i];

		return (float) sum;
	}

	// finds the maximum element in a vector
	static short findMax(float* vec, int n)
	{
		float maxVec = 0;
		short idx= 0;
		short i;

		for (i = 1; i < n - 1; ++i) {
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
			hist[c.GetR()]++;
			hist[c.GetG()]++;
			hist[c.GetB()]++;
			if (hasSemiTransparency)
				hist[c.GetA()]++;
		}
	}
		
	short getOtsuThreshold(const vector<ARGB>& pixels)
	{
	    float vet[256] = { 0 };
		int hist[256] = { 0 };
		
		getHistogram(pixels, hist);

		// loop through all possible t values and maximize between class variance
		for (int k = 1; k != BYTE_MAX; k++) {
			float p1 = Px(0, k, hist);
			float p2 = Px(k + 1, BYTE_MAX, hist);
			float p12 = p1 * p2;
			if (p12 == 0) 
				p12 = 1;
			float diff = (Mx(0, k, hist) * p2) - (Mx(k + 1, BYTE_MAX, hist) * p1);
			vet[k] = (float)diff * diff / p12;
		}

		return findMax(vet, 256);
	}
	
	bool threshold(vector<ARGB>& pixels, short thresh, float weight = 1.0f)
	{
		if (thresh >= 200)
		{
			weight = .75f;
			thresh = 200;
		}

		auto minThresh = (BYTE)(thresh * weight);
		auto maxThresh = (BYTE)thresh;
		for (int i = 0; i < pixels.size(); ++i) {
			Color c(pixels[i]);
			if (c.GetR() + c.GetG() + c.GetB() > maxThresh * 3)
				pixels[i] = Color::MakeARGB(c.GetA(), BYTE_MAX, BYTE_MAX, BYTE_MAX);
			else if (c.GetR() + c.GetG() + c.GetB() < minThresh * 3)
				pixels[i] = Color::MakeARGB(c.GetA(), 0, 0, 0);
		}
		
		return true;
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

		double mindist = INT_MAX;
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
				auto grey = min(ptr[0], ptr[1]);
				grey = min(ptr[1], ptr[2]);
				if (min1 > grey)
					min1 = grey;

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
				auto grey = min(ptr[0], ptr[1]);
				grey = min(ptr[1], ptr[2]);
				ptr[0] = ptr[1] = ptr[2] = (BYTE)((grey - min1) * (BYTE_MAX / (max1 - min1)));
				ptr += DJ;
			}
			ptr += remain;
		}

		pSourceImg->UnlockBits(&data);
		return pSourceImg;
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	bool Otsu::ConvertGrayScaleToBinary(Bitmap* pSrcImg, Bitmap* pDest)
	{
		auto pSourceImg = unique_ptr<Bitmap>(ConvertToGrayScale(pSrcImg));

		auto bitmapWidth = pSourceImg->GetWidth();
		auto bitmapHeight = pSourceImg->GetHeight();

		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		if (!GrabPixels(pSourceImg.get(), pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor))
			return false;

		auto otsuThreshold = getOtsuThreshold(pixels);
		if (!threshold(pixels, otsuThreshold))
			return false;

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