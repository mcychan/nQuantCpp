#pragma once
/**
 * A blue-noise-based dither does not diffuse error, and uses a tiling blue noise pattern, 
 * with a fine-grained checker board pattern
 * and a roughly-white-noise pattern obtained by distorting the blue noise, but only applies these noisy pattern
 * when there's error matching a color from the image to a color in the palette. 
 * Copyright (c) 2021 Miller Cy Chan */

#include "stdafx.h"
#include "BlueNoise.h"

#include <memory>

namespace BlueNoise
{	
	
	void dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, const float weight)
	{
		const float strength = 1 / 3.0f;  	
		
		for (UINT y = 0; y < height; ++y) {
			for (UINT x = 0; x < width; ++x) {
				Color pixel(pixels[x + y * width]);
				int r_pix = pixel.GetR();
				int g_pix = pixel.GetG();
				int b_pix = pixel.GetB();
				int a_pix = pixel.GetA();

				Color c1 = pPalette->Entries[qPixels[x + y * width]];
				float adj = (RAW_BLUE_NOISE[(x & 63) | (y & 63) << 6] + 0.5f) / 127.5f;
				adj += ((x + y & 1) - 0.5f) * strength / 8.0f;
				adj *= weight;
				r_pix = static_cast<BYTE>(min(BYTE_MAX, max(r_pix + (adj * (r_pix - c1.GetR())), 0)));
				g_pix = static_cast<BYTE>(min(BYTE_MAX, max(g_pix + (adj * (g_pix - c1.GetG())), 0)));
				b_pix = static_cast<BYTE>(min(BYTE_MAX, max(b_pix + (adj * (b_pix - c1.GetB())), 0)));
				a_pix = static_cast<BYTE>(min(BYTE_MAX, max(a_pix + (adj * (a_pix - c1.GetA())), 0)));

				c1 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
				qPixels[x + y * width] = ditherFn(pPalette, c1.GetValue(), x + y);
			}
		}
    }
}
