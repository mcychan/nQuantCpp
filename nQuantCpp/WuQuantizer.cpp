#pragma once
///////////////////////////////////////////////////////////////////////
//	    C Implementation of Wu's Color Quantizer (v. 2)
//	    (see Graphics Gems vol. II, pp. 126-133)
//
// Author:	Xiaolin Wu
// Dept. of Computer Science
// Univ. of Western Ontario
// London, Ontario N6A 5B7
// wu@csd.uwo.ca
// 
// Algorithm: Greedy orthogonal bipartition of RGB space for variance
// 	   minimization aided by inclusion-exclusion tricks.
// 	   For speed no nearest neighbor search is done. Slightly
// 	   better performance can be expected by more sophisticated
// 	   but more expensive versions.
// 
// Free to distribute, comments and suggestions are appreciated.
///////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "WuQuantizer.h"
#include <map>

namespace nQuant
{
	/// <summary><para>Shift color values right this many bits.</para><para>This reduces the granularity of the color maps produced, making it much faster.</para></summary>
	/// 3 = value error of 8 (0 and 7 will look the same to it, 0 and 8 different); Takes ~4MB for color tables; ~.25 -> .50 seconds
	/// 2 = value error of 4; Takes ~64MB for color tables; ~3 seconds
	/// RAM usage roughly estimated with: ( ( 256 >> SidePixShift ) ^ 4 ) * 60
	/// Default SidePixShift = 3
	const byte SIDEPIXSHIFT = 3;
	const byte MAXSIDEINDEX = 256 / (1 << SIDEPIXSHIFT);
	const byte SIDESIZE = MAXSIDEINDEX + 1;
	const UINT TOTAL_SIDESIZE = SIDESIZE * SIDESIZE * SIDESIZE * SIDESIZE;
	
	bool hasTransparency = false;
	ARGB m_transparentColor;
	map<ARGB, vector<short> > closestMap;
	map<ARGB, UINT> rightMatches;

	struct Box {
		byte AlphaMinimum = 0;
		byte AlphaMaximum = 0;
		byte RedMinimum = 0;
		byte RedMaximum = 0;
		byte GreenMinimum = 0;
		byte GreenMaximum = 0;
		byte BlueMinimum = 0;
		byte BlueMaximum = 0;
		long Size = 0;
	};

	struct LookupData {
		vector<ARGB> lookups;
		unique_ptr<UINT[]> tags;

		LookupData(UINT sideSize) {
			tags = make_unique<UINT[]>(TOTAL_SIDESIZE);
		}
	};

	struct QuantizedPalette {
	private:		
		unique_ptr<byte[]> pPaletteBytes;
		ColorPalette* pPalette = nullptr;

	public:
		unique_ptr<UINT[]> pixelIndex;
		QuantizedPalette(UINT size, const int nMaxColors) {
			pixelIndex = make_unique<UINT[]>(size);
			pPaletteBytes = make_unique<byte[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
			pPalette = (ColorPalette*) pPaletteBytes.get();
			pPalette->Count = nMaxColors;
		}

		inline ColorPalette* GetPalette() const {
			return pPalette;
		}
	};

	struct CubeCut {
		bool valid;
		BYTE position;
		double value;

		CubeCut(bool isValid, BYTE cutPoint, double result) {
			valid = isValid;
			position = cutPoint;
			value = result;
		}
	};

	struct ColorData {
		unique_ptr<long[]> weights;
		unique_ptr<long[]> momentsAlpha;
		unique_ptr<long[]> momentsRed;
		unique_ptr<long[]> momentsGreen;
		unique_ptr<long[]> momentsBlue;
		unique_ptr<double[]> moments;

		unique_ptr<ARGB[]> pixels;
		unique_ptr<ARGB[]> quantizedPixels;

		UINT pixelsCount = 0;
		UINT pixelFillingCounter = 0;

		ColorData(UINT sideSize, UINT bitmapWidth, UINT bitmapHeight) {
			const int TOTAL_SIDESIZE = sideSize * sideSize * sideSize * sideSize;
			weights = make_unique<long[]>(TOTAL_SIDESIZE);
			momentsAlpha = make_unique<long[]>(TOTAL_SIDESIZE);
			momentsRed = make_unique<long[]>(TOTAL_SIDESIZE);
			momentsGreen = make_unique<long[]>(TOTAL_SIDESIZE);
			momentsBlue = make_unique<long[]>(TOTAL_SIDESIZE);
			moments = make_unique<double[]>(TOTAL_SIDESIZE);
			pixelsCount = bitmapWidth * bitmapHeight;
			pixels = make_unique<ARGB[]>(pixelsCount);
			quantizedPixels = make_unique<ARGB[]>(pixelsCount);
		}
		
		inline ARGB* GetPixels() {
			return pixels.get();
		}

		void AddPixel(ARGB pixel, int quantizedPixel)
		{
			pixels[pixelFillingCounter] = pixel;
			quantizedPixels[pixelFillingCounter++] = quantizedPixel;
		}
	};

	inline UINT Index(byte alpha, byte red, byte green, byte blue) {
		return alpha + red * SIDESIZE + green * SIDESIZE * SIDESIZE + blue * SIDESIZE * SIDESIZE * SIDESIZE;
	}

	inline long Volume(const Box& cube, long* moment)
	{
		return (moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
			- (moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);
	}

	inline double Volume(const Box& cube, double* moment)
	{
		return (moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
			- (moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);
	}

	inline long Top(const Box& cube, Pixel direction, byte position, long* moment)
	{
		switch (direction)
		{
		case Alpha:
			return (moment[Index(position, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(position, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(position, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(position, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
				- (moment[Index(position, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(position, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(position, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(position, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);

		case Red:
			return (moment[Index(cube.AlphaMaximum, position, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, position, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, position, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, position, cube.GreenMinimum, cube.BlueMaximum)])
				- (moment[Index(cube.AlphaMaximum, position, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, position, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, position, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, position, cube.GreenMinimum, cube.BlueMinimum)]);

		case Green:
			return (moment[Index(cube.AlphaMaximum, cube.RedMaximum, position, cube.BlueMaximum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, position, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, position, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, position, cube.BlueMaximum)])
				- (moment[Index(cube.AlphaMaximum, cube.RedMaximum, position, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, position, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, position, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, position, cube.BlueMinimum)]);

		case Blue:
			return (moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, position)] - moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, position)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, position)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, position)])
				- (moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, position)] - moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, position)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, position)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, position)]);

		default:
			return 0;
		}
	}

	inline long Bottom(const Box& cube, Pixel direction, long* moment)
	{
		switch (direction)
		{
		case Alpha:
			return -(moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
				- (-moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);

		case Red:
			return -(moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
				- (-moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);

		case Green:
			return -(moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMaximum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMaximum)])
				- (-moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);

		case Blue:
			return -(moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMaximum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)])
				- (-moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMaximum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMaximum, cube.GreenMinimum, cube.BlueMinimum)] + moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMaximum, cube.BlueMinimum)] - moment[Index(cube.AlphaMinimum, cube.RedMinimum, cube.GreenMinimum, cube.BlueMinimum)]);

		default:
			return 0;
		}
	}	

	void CompileColorData(ColorData& colorData, const Color& color, const byte alphaThreshold, const byte alphaFader)
	{
		byte pixelBlue = color.GetBlue();
		byte pixelGreen = color.GetGreen();
		byte pixelRed = color.GetRed();
		byte pixelAlpha = color.GetAlpha();

		byte indexAlpha = (pixelAlpha >> SIDEPIXSHIFT) + 1;
		byte indexRed = (pixelRed >> SIDEPIXSHIFT) + 1;
		byte indexGreen = (pixelGreen >> SIDEPIXSHIFT) + 1;
		byte indexBlue = (pixelBlue >> SIDEPIXSHIFT) + 1;

		if (pixelAlpha > alphaThreshold) {
			if (pixelAlpha < BYTE_MAX) {
				short alpha = pixelAlpha + (pixelAlpha % alphaFader);
				pixelAlpha = alpha > BYTE_MAX ? BYTE_MAX : alpha;
				indexAlpha = (pixelAlpha >> 3) + 1;
			}

			const int index = Index(indexAlpha, indexRed, indexGreen, indexBlue);
			if (index < TOTAL_SIDESIZE) {
				colorData.weights[index]++;
				colorData.momentsRed[index] += pixelRed;
				colorData.momentsGreen[index] += pixelGreen;
				colorData.momentsBlue[index] += pixelBlue;
				colorData.momentsAlpha[index] += pixelAlpha;
				colorData.moments[index] += (pixelAlpha * pixelAlpha) + (pixelRed * pixelRed) + (pixelGreen * pixelGreen) + (pixelBlue * pixelBlue);
			}
		}

		ARGB pixel = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
		ARGB qPixel = Color::MakeARGB(indexAlpha, indexRed, indexGreen, indexBlue);
		colorData.AddPixel(pixel, qPixel);
	}

	void BuildHistogram(ColorData& colorData, Bitmap* sourceImage, byte alphaThreshold, byte alphaFader)
	{
		UINT bitDepth = GetPixelFormatSize(sourceImage->GetPixelFormat());
		UINT bitmapWidth = sourceImage->GetWidth();
		UINT bitmapHeight = sourceImage->GetHeight();

		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++) {
					Color color;
					sourceImage->GetPixel(x, y, &color);
					if(color.GetA() < BYTE_MAX)
						hasTransparency = true;
					if(color.GetA() == 0)
						m_transparentColor = color.GetValue();
					CompileColorData(colorData, color, alphaThreshold, alphaFader);
				}
			}
			return;
		}

		BitmapData data;
		Status status = sourceImage->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, sourceImage->GetPixelFormat(), &data);
		if (status != Ok)
			return;		

		auto pRowSource = (byte*) data.Scan0;
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

		// First loop: gather color information
		for (UINT y = 0; y < bitmapHeight; y++)	// For each row...
		{
			auto pPixelSource = pRowSource;

			for (UINT x = 0; x < bitmapWidth; x++)	// ...for each pixel...
			{
				byte pixelBlue = *pPixelSource++;
				byte pixelGreen = *pPixelSource++;
				byte pixelRed = *pPixelSource++;
				byte pixelAlpha = bitDepth < 32 ? BYTE_MAX : *pPixelSource++;
				if (pixelAlpha < BYTE_MAX)
					hasTransparency = true;
				
				Color color(Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue));
				if(pixelAlpha == 0)
					m_transparentColor = color.GetValue();
				CompileColorData(colorData, color, alphaThreshold, alphaFader);
			}

			pRowSource += strideSource;
		}

		sourceImage->UnlockBits(&data);
	}

	void CalculateMoments(ColorData& data)
	{
		for (int alphaIndex = 1; alphaIndex <= MAXSIDEINDEX; ++alphaIndex)
		{
			UINT xarea[SIDESIZE][SIDESIZE][SIDESIZE] = { 0 };
			UINT xareaAlpha[SIDESIZE][SIDESIZE][SIDESIZE] = { 0 };
			UINT xareaRed[SIDESIZE][SIDESIZE][SIDESIZE] = { 0 };
			UINT xareaGreen[SIDESIZE][SIDESIZE][SIDESIZE] = { 0 };
			UINT xareaBlue[SIDESIZE][SIDESIZE][SIDESIZE] = { 0 };
			double xarea2[SIDESIZE][SIDESIZE][SIDESIZE] = { 0.0 };
						
			for (int redIndex = 1; redIndex <= MAXSIDEINDEX; ++redIndex)
			{
				UINT area[SIDESIZE] = { 0 };
				UINT areaAlpha[SIDESIZE] = { 0 };
				UINT areaRed[SIDESIZE] = { 0 };
				UINT areaGreen[SIDESIZE] = { 0 };
				UINT areaBlue[SIDESIZE] = { 0 };
				double area2[SIDESIZE] = { 0.0 };

				for (int greenIndex = 1; greenIndex <= MAXSIDEINDEX; ++greenIndex) {
					volatile UINT line = 0;
					volatile UINT lineAlpha = 0;
					volatile UINT lineRed = 0;
					volatile UINT lineGreen = 0;
					volatile UINT lineBlue = 0;
					volatile double line2 = 0.0;

					for (int blueIndex = 1; blueIndex <= MAXSIDEINDEX; ++blueIndex) {
						line += data.weights[Index(alphaIndex, redIndex, greenIndex, blueIndex)];
						lineAlpha += data.momentsAlpha[Index(alphaIndex, redIndex, greenIndex, blueIndex)];
						lineRed += data.momentsRed[Index(alphaIndex, redIndex, greenIndex, blueIndex)];
						lineGreen += data.momentsGreen[Index(alphaIndex, redIndex, greenIndex, blueIndex)];
						lineBlue += data.momentsBlue[Index(alphaIndex, redIndex, greenIndex, blueIndex)];
						line2 += data.moments[Index(alphaIndex, redIndex, greenIndex, blueIndex)];

						area[blueIndex] += line;
						areaAlpha[blueIndex] += lineAlpha;
						areaRed[blueIndex] += lineRed;
						areaGreen[blueIndex] += lineGreen;
						areaBlue[blueIndex] += lineBlue;
						area2[blueIndex] += line2;

						xarea[redIndex][greenIndex][blueIndex] = xarea[redIndex - 1][greenIndex][blueIndex] + area[blueIndex];
						xareaAlpha[redIndex][greenIndex][blueIndex] = xareaAlpha[redIndex - 1][greenIndex][blueIndex] + areaAlpha[blueIndex];
						xareaRed[redIndex][greenIndex][blueIndex] = xareaRed[redIndex - 1][greenIndex][blueIndex] + areaRed[blueIndex];
						xareaGreen[redIndex][greenIndex][blueIndex] = xareaGreen[redIndex - 1][greenIndex][blueIndex] + areaGreen[blueIndex];
						xareaBlue[redIndex][greenIndex][blueIndex] = xareaBlue[redIndex - 1][greenIndex][blueIndex] + areaBlue[blueIndex];
						xarea2[redIndex][greenIndex][blueIndex] = xarea2[redIndex - 1][greenIndex][blueIndex] + area2[blueIndex];

						data.weights[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.weights[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xarea[redIndex][greenIndex][blueIndex];
						data.momentsAlpha[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.momentsAlpha[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xareaAlpha[redIndex][greenIndex][blueIndex];
						data.momentsRed[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.momentsRed[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xareaRed[redIndex][greenIndex][blueIndex];
						data.momentsGreen[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.momentsGreen[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xareaGreen[redIndex][greenIndex][blueIndex];
						data.momentsBlue[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.momentsBlue[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xareaBlue[redIndex][greenIndex][blueIndex];
						data.moments[Index(alphaIndex, redIndex, greenIndex, blueIndex)] = data.moments[Index(alphaIndex - 1, redIndex, greenIndex, blueIndex)] + xarea2[redIndex][greenIndex][blueIndex];
					}
				}
			}
		}
	}	

	CubeCut Maximize(ColorData& data, const Box& cube, Pixel direction, byte first, byte last, long wholeAlpha, long wholeRed, long wholeGreen, long wholeBlue, long wholeWeight)
	{
		auto bottomAlpha = Bottom(cube, direction, data.momentsAlpha.get());
		auto bottomRed = Bottom(cube, direction, data.momentsRed.get());
		auto bottomGreen = Bottom(cube, direction, data.momentsGreen.get());
		auto bottomBlue = Bottom(cube, direction, data.momentsBlue.get());
		auto bottomWeight = Bottom(cube, direction, data.weights.get());

		volatile bool valid = false;
		volatile double result = 0.0;
		volatile byte cutPoint = 0;

		#pragma omp parallel for
		for (int position = first; position < last; ++position)
		{
			long halfAlpha = bottomAlpha + Top(cube, direction, position, data.momentsAlpha.get());
			long halfRed = bottomRed + Top(cube, direction, position, data.momentsRed.get());
			long halfGreen = bottomGreen + Top(cube, direction, position, data.momentsGreen.get());
			long halfBlue = bottomBlue + Top(cube, direction, position, data.momentsBlue.get());
			long halfWeight = bottomWeight + Top(cube, direction, position, data.weights.get());

			if (halfWeight == 0)
				continue;

			long halfDistance = halfAlpha * halfAlpha + halfRed * halfRed + halfGreen * halfGreen + halfBlue * halfBlue;
			double temp = halfDistance * 1.0 / halfWeight;

			halfAlpha = wholeAlpha - halfAlpha;
			halfRed = wholeRed - halfRed;
			halfGreen = wholeGreen - halfGreen;
			halfBlue = wholeBlue - halfBlue;
			halfWeight = wholeWeight - halfWeight;

			if (halfWeight != 0) {
				halfDistance = halfAlpha * halfAlpha + halfRed * halfRed + halfGreen * halfGreen + halfBlue * halfBlue;
				temp += halfDistance / halfWeight;

				if (temp > result) {
					valid = true;
					result = temp;
					cutPoint = position;
				}
			}
		}

		return CubeCut(valid, cutPoint, result);
	}

	bool Cut(ColorData& data, Box& first, Box& second)
	{
		Pixel direction = Blue;
		auto wholeAlpha = Volume(first, data.momentsAlpha.get());
		auto wholeRed = Volume(first, data.momentsRed.get());
		auto wholeGreen = Volume(first, data.momentsGreen.get());
		auto wholeBlue = Volume(first, data.momentsBlue.get());
		auto wholeWeight = Volume(first, data.weights.get());

		auto maxAlpha = Maximize(data, first, Alpha, static_cast<byte>(first.AlphaMinimum + 1), first.AlphaMaximum, wholeAlpha, wholeRed, wholeGreen, wholeBlue, wholeWeight);
		auto maxRed = Maximize(data, first, Red, static_cast<byte>(first.RedMinimum + 1), first.RedMaximum, wholeAlpha, wholeRed, wholeGreen, wholeBlue, wholeWeight);
		auto maxGreen = Maximize(data, first, Green, static_cast<byte>(first.GreenMinimum + 1), first.GreenMaximum, wholeAlpha, wholeRed, wholeGreen, wholeBlue, wholeWeight);
		auto maxBlue = Maximize(data, first, Blue, static_cast<byte>(first.BlueMinimum + 1), first.BlueMaximum, wholeAlpha, wholeRed, wholeGreen, wholeBlue, wholeWeight);

		second.AlphaMaximum = first.AlphaMaximum;
		second.RedMaximum = first.RedMaximum;
		second.GreenMaximum = first.GreenMaximum;
		second.BlueMaximum = first.BlueMaximum;

		if ((maxAlpha.value >= maxRed.value) && (maxAlpha.value >= maxGreen.value) && (maxAlpha.value >= maxBlue.value)) {
			if (!maxAlpha.valid)
				return false;
			direction = Alpha;
		}
		else if ((maxRed.value >= maxAlpha.value) && (maxRed.value >= maxGreen.value) && (maxRed.value >= maxBlue.value))
			direction = Red;
		else if ((maxGreen.value >= maxAlpha.value) && (maxGreen.value >= maxRed.value) && (maxGreen.value >= maxBlue.value))
			direction = Green;

		switch (direction)
		{
		case Alpha:
			second.AlphaMinimum = first.AlphaMaximum = maxAlpha.position;
			second.RedMinimum = first.RedMinimum;
			second.GreenMinimum = first.GreenMinimum;
			second.BlueMinimum = first.BlueMinimum;
			break;

		case Red:
			second.RedMinimum = first.RedMaximum = maxRed.position;
			second.AlphaMinimum = first.AlphaMinimum;
			second.GreenMinimum = first.GreenMinimum;
			second.BlueMinimum = first.BlueMinimum;
			break;

		case Green:
			second.GreenMinimum = first.GreenMaximum = maxGreen.position;
			second.AlphaMinimum = first.AlphaMinimum;
			second.RedMinimum = first.RedMinimum;
			second.BlueMinimum = first.BlueMinimum;
			break;

		case Blue:
			second.BlueMinimum = first.BlueMaximum = maxBlue.position;
			second.AlphaMinimum = first.AlphaMinimum;
			second.RedMinimum = first.RedMinimum;
			second.GreenMinimum = first.GreenMinimum;
			break;
		}

		first.Size = (first.AlphaMaximum - first.AlphaMinimum) * (first.RedMaximum - first.RedMinimum) * (first.GreenMaximum - first.GreenMinimum) * (first.BlueMaximum - first.BlueMinimum);
		second.Size = (second.AlphaMaximum - second.AlphaMinimum) * (second.RedMaximum - second.RedMinimum) * (second.GreenMaximum - second.GreenMinimum) * (second.BlueMaximum - second.BlueMinimum);

		return true;
	}

	double CalculateVariance(ColorData& data, const Box& cube)
	{
		auto volumeAlpha = Volume(cube, data.momentsAlpha.get());
		auto volumeRed = Volume(cube, data.momentsRed.get());
		auto volumeGreen = Volume(cube, data.momentsGreen.get());
		auto volumeBlue = Volume(cube, data.momentsBlue.get());
		auto volumeMoment = Volume(cube, data.moments.get());
		auto volumeWeight = Volume(cube, data.weights.get());

		long distance = volumeAlpha * volumeAlpha + volumeRed * volumeRed + volumeGreen * volumeGreen + volumeBlue * volumeBlue;

		return volumeWeight != 0 ? (volumeMoment - distance / volumeWeight) : 0.0;
	}	

	void SplitData(vector<Box>& boxList, UINT& colorCount, ColorData& data)
	{
		const int COLORSIZE = --colorCount;
		int next = 0;
		auto volumeVariance = make_unique<double[]>(COLORSIZE);
		boxList.resize(COLORSIZE);
		boxList[0].AlphaMaximum = MAXSIDEINDEX;
		boxList[0].RedMaximum = MAXSIDEINDEX;
		boxList[0].GreenMaximum = MAXSIDEINDEX;
		boxList[0].BlueMaximum = MAXSIDEINDEX;

		for (int cubeIndex = 1; cubeIndex < COLORSIZE; ++cubeIndex) {
			if (Cut(data, boxList[next], boxList[cubeIndex])) {
				volumeVariance[next] = boxList[next].Size > 1 ? CalculateVariance(data, boxList[next]) : 0.0;
				volumeVariance[cubeIndex] = boxList[cubeIndex].Size > 1 ? CalculateVariance(data, boxList[cubeIndex]) : 0.0;
			}
			else {
				volumeVariance[next] = 0.0;
				cubeIndex--;
			}

			next = 0;
			double temp = volumeVariance[0];

			for (int index = 1; index <= cubeIndex; ++index) {
				if (volumeVariance[index] <= temp)
					continue;
				temp = volumeVariance[index];
				next = index;
			}

			if (temp > 0.0)
				continue;

			colorCount = cubeIndex + 1;
			break;
		}
	}

	void BuildLookups(LookupData& lookupData, vector<Box>& cubes, const ColorData& data)
	{
		volatile UINT lookupsCount = 0;
		for (auto const& cube : cubes) {
			#pragma omp parallel for
			for (int alphaIndex = cube.AlphaMinimum + 1; alphaIndex <= cube.AlphaMaximum; ++alphaIndex) {
				for (int redIndex = cube.RedMinimum + 1; redIndex <= cube.RedMaximum; ++redIndex) {
					for (int greenIndex = cube.GreenMinimum + 1; greenIndex <= cube.GreenMaximum; ++greenIndex) {
						for (int blueIndex = cube.BlueMinimum + 1; blueIndex <= cube.BlueMaximum; ++blueIndex) {
							const auto index = Index(alphaIndex, redIndex, greenIndex, blueIndex);
							if (index < TOTAL_SIDESIZE)
								lookupData.tags[index] = lookupsCount;
						}
					}
				}
			}

			auto weight = Volume(cube, data.weights.get());

			if (weight <= 0)
				continue;

			byte alpha = static_cast<byte>(Volume(cube, data.momentsAlpha.get()) / weight);
			byte red = static_cast<byte>(Volume(cube, data.momentsRed.get()) / weight);
			byte green = static_cast<byte>(Volume(cube, data.momentsGreen.get()) / weight);
			byte blue = static_cast<byte>(Volume(cube, data.momentsBlue.get()) / weight);
			lookupData.lookups.emplace_back(Color::MakeARGB(alpha, red, green, blue));
			++lookupsCount;
		}
	}
	
	UINT bestcolor(const LookupData& lookupData, ARGB argb, int pixelIndex, byte alphaThreshold)
	{
		Color c(argb);
		UINT k = 0;
		if (c.GetA() <= alphaThreshold)			
			return k;

		if (hasTransparency) {
			auto got = rightMatches.find(argb);
			if (got == rightMatches.end()) {
				int bestDistance = INT_MAX;
				auto lookups = lookupData.lookups;
				auto lookupsCount = lookups.size();

				for (size_t lookupIndex = 0; lookupIndex < lookupsCount; lookupIndex++) {
					Color lookup(lookups[lookupIndex]);
					int deltaAlpha = abs(c.GetA() - lookup.GetA());
					int deltaRed = abs(c.GetR() - lookup.GetR());
					int deltaGreen = abs(c.GetG() - lookup.GetG());
					int deltaBlue = abs(c.GetB() - lookup.GetB());

					int distance = deltaAlpha + deltaRed + deltaGreen + deltaBlue;

					if (distance >= bestDistance)
						continue;

					bestDistance = distance;
					k = lookupIndex;
				}

				rightMatches[argb] = k;
			}
			else
				k = got->second;

			return k;
		}

		vector<short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHORT_MAX;

			auto lookups = lookupData.lookups;
			auto lookupsCount = lookups.size();

			for (; k < lookupsCount; k++) {
				Color c2(lookups[k]);
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

			if (closest[3] == 100000000)
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
	
	void GetQuantizedPalette(QuantizedPalette& qPalette, const ColorData& data, const LookupData& lookupData, UINT colorCount, byte alphaThreshold)
	{
		const UINT COLOR_SIZE = colorCount + 1;
		auto alphas = make_unique<UINT[]>(COLOR_SIZE);
		auto reds = make_unique<UINT[]>(COLOR_SIZE);
		auto greens = make_unique<UINT[]>(COLOR_SIZE);
		auto blues = make_unique<UINT[]>(COLOR_SIZE);
		auto sums = make_unique<UINT[]>(COLOR_SIZE);

		int pixelsCount = data.pixelsCount;

		for (int pixelIndex = 0; pixelIndex < pixelsCount; pixelIndex++) {
			Color color(data.quantizedPixels[pixelIndex]);
			data.quantizedPixels[pixelIndex] = lookupData.tags[Index(color.GetA(), color.GetR(), color.GetG(), color.GetB())];

			auto argb = data.pixels[pixelIndex];			

			UINT bestMatch = bestcolor(lookupData, argb, pixelIndex, alphaThreshold);

			if (bestMatch < COLOR_SIZE) {
				Color pixel(argb);
				alphas[bestMatch] += pixel.GetA();
				reds[bestMatch] += pixel.GetR();
				greens[bestMatch] += pixel.GetG();
				blues[bestMatch] += pixel.GetB();
				sums[bestMatch]++;
			}

			qPalette.pixelIndex[pixelIndex] = bestMatch;
		}

		UINT paletteIndex = 0;
		if (hasTransparency) {
			++paletteIndex;
			--colorCount;
		}

		for (; paletteIndex < colorCount; paletteIndex++) {
			if (sums[paletteIndex] <= 0)
				continue;

			alphas[paletteIndex] /= sums[paletteIndex];
			reds[paletteIndex] /= sums[paletteIndex];
			greens[paletteIndex] /= sums[paletteIndex];
			blues[paletteIndex] /= sums[paletteIndex];

			auto color = Color::MakeARGB(alphas[paletteIndex], reds[paletteIndex], greens[paletteIndex], blues[paletteIndex]);
			qPalette.GetPalette()->Entries[paletteIndex] = color;
		}

		if(hasTransparency)
			qPalette.GetPalette()->Entries[0] = m_transparentColor;
	}

	bool quantize_image(const ARGB* pixels, const LookupData& lookupData, UINT* qPixels, UINT nMaxColors, int width, int height, bool dither, byte alphaThreshold)
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
					qPixels[pixelIndex] = bestcolor(lookupData, argb, i, alphaThreshold);

					auto lookups = lookupData.lookups;
					Color c2(lookups[qPixels[pixelIndex]]);
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
				qPixels[i] = bestcolor(lookupData, pixels[i], i, alphaThreshold);
		}
		return true;
	}

	bool ProcessImagePixels(Bitmap* pDest, const QuantizedPalette& qPalette)
	{
		auto palette = qPalette.GetPalette();
		pDest->SetPalette(palette);
		UINT nMaxColors = palette->Count;

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
				byte index = static_cast<byte>(qPalette.pixelIndex[pixelIndex++]);

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
				}
			}

			pRowDest += strideDest;
		}

		status = pDest->UnlockBits(&targetData);
		return pDest->GetLastStatus() == Ok;
	}

	bool WuQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither, byte alphaThreshold, byte alphaFader)
	{
		hasTransparency = false;
		ColorData colorData(SIDESIZE, pSource->GetWidth(), pSource->GetHeight());
		BuildHistogram(colorData, pSource, alphaThreshold, alphaFader);
		CalculateMoments(colorData);
		vector<Box> cubes;
		SplitData(cubes, nMaxColors, colorData);		
		LookupData lookupData(SIDESIZE);
		BuildLookups(lookupData, cubes, colorData);
		cubes.clear();
		QuantizedPalette qPalette(colorData.pixelsCount, nMaxColors);
		GetQuantizedPalette(qPalette, colorData, lookupData, nMaxColors, alphaThreshold);
		
		quantize_image(colorData.GetPixels(), lookupData, qPalette.pixelIndex.get(), nMaxColors, pSource->GetWidth(), pSource->GetHeight(), dither, alphaThreshold);
		closestMap.clear();
		rightMatches.clear();

		return ProcessImagePixels(pDest, qPalette);
	}

}
