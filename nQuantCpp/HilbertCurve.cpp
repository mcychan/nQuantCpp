#pragma once
/* The Hilbert curve is a space filling curve that visits every point in a square grid with a size of any other power of 2.
Copyright (c) 2021 Miller Cy Chan
* It was first described by David Hilbert in 1892. Applications of the Hilbert curve are in image processing: especially image compression and dithering. */

#include "stdafx.h"
#include "HilbertCurve.h"

#include <memory>
#include <vector>

namespace Riemersma
{
	enum Direction { LEFT, RIGHT, DOWN, UP };
	
    struct ErrorBox
    {
        float p[4] = { 0 };

        ErrorBox() {
        }
        ErrorBox(const Color& c) {
            p[0] = c.GetR();
            p[1] = c.GetG();
            p[2] = c.GetB();
            p[3] = c.GetA();
        }
        inline float& operator[](int index)
        {
            return p[index];
        }
    };
	
	int x, y;
	UINT m_width, m_height;
	const ARGB* m_image;
	const ColorPalette* m_pPalette;
	unsigned short* m_qPixels;
	DitherFn m_ditherFn;
	vector<ErrorBox> errorq;
	float* m_weights;
    
	static const BYTE DITHER_MAX = 16;
	static const float BLOCK_SIZE = 256.0f; 
    
    void ditherCurrentPixel()
	{
	    if(x >= 0 && y >= 0 && x < m_width && y < m_height) {
	    	Color pixel(m_image[x + y * m_width]);
	    	ErrorBox error(pixel);	    	
	        for(int c = 0; c < DITHER_MAX; ++c) {
	        	auto& eb = errorq[c];
	        	for(int j = 0; j < sizeof(eb.p) / sizeof(float); ++j)
	        		error[j] += eb[j] * m_weights[c];
	        }

            auto r_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[0], 0)));
            auto g_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[1], 0)));
            auto b_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[2], 0)));
            auto a_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[3], 0)));
	        
	        Color c1 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);		        
            m_qPixels[x + y * m_width] = m_ditherFn(m_pPalette, m_pPalette->Count, c1.GetValue());

	        errorq.erase(errorq.begin());
	        Color c2 = m_pPalette->Entries[m_qPixels[x + y * m_width]];
	        error[0] = r_pix - c2.GetR();
	        error[1] = g_pix - c2.GetG();
	        error[2] = b_pix - c2.GetB();
	        error[3] = a_pix - c2.GetA();
	        
	        for(int j = 0; j < sizeof(error.p) / sizeof(float); ++j) {
	        	if(abs(error[j]) > DITHER_MAX)
	        		error[j] = error[j] < 0 ? DITHER_MAX : DITHER_MAX;
	        }
	        errorq.emplace_back(error);
	    }
	}
    
    void navTo(Direction dir)
	{
    	ditherCurrentPixel();
		switch(dir)
        {
            case LEFT:
            	--x;
            	break;
            case RIGHT:
            	++x;
            	break;
            case UP:
            	--y;
            	break;
            case DOWN:
            	++y;
            	break;
        }
	}

    void curve(const int level, Direction a, Direction b, Direction c, Direction d, Direction e, Direction f, Direction g)
    {
        auto iter = [](const int level, Direction dir)
        {
            if (level <= 0)
                return;

            switch (dir)
            {
            case LEFT:
                curve(level, UP, LEFT, LEFT, DOWN, RIGHT, DOWN, LEFT);
                break;
            case RIGHT:
                curve(level, DOWN, RIGHT, RIGHT, UP, LEFT, UP, RIGHT);
                break;
            case UP:
                curve(level, LEFT, UP, UP, RIGHT, DOWN, RIGHT, UP);
                break;
            case DOWN:
                curve(level, RIGHT, DOWN, DOWN, LEFT, UP, LEFT, DOWN);
                break;
            }
        };

		iter(level-1, a);
		navTo(e);
		iter(level-1, b);
		navTo(f);
        iter(level-1, c);
        navTo(g);
        iter(level-1, d);
    }    
	
	void HilbertCurve::dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, unsigned short* qPixels)
    {
		m_width = width;
        m_height = height;
        m_image = pixels;
        m_pPalette = pPalette;
        m_qPixels = qPixels;
        m_ditherFn = ditherFn;
        auto pWeights = make_unique<float[]>(DITHER_MAX);
        m_weights = pWeights.get();
    	
        /* Dithers all pixels of the image in sequence using
         * the Hilbert path, and distributes the error in
         * a sequence of 16 pixels.
         */
        x = y = 0;
        errorq.clear();
        const float weightRatio = (float)pow(BLOCK_SIZE + 1.0f, 1.0f / (DITHER_MAX - 1.0f));
        float weight = 1.0f, sumweight = 0.0f;
        for (int c = 0; c < DITHER_MAX; ++c)
        {
            errorq.resize(DITHER_MAX);
            sumweight += (m_weights[DITHER_MAX - c - 1] = 1.0f / weight);
            weight *= weightRatio;
        }

        weight = 0.0f; /* Normalize */
        for (int c = 0; c < DITHER_MAX; ++c)
            weight += (m_weights[c] /= sumweight);
        m_weights[0] += 1.0f - weight;
        /* Walk the path. */
        int i = max(m_width, m_height), depth = 0;
        while (i > 0) {
            ++depth;
            i >>= 1;
        }

        curve(depth, UP, LEFT, LEFT, DOWN, RIGHT, DOWN, LEFT);
        ditherCurrentPixel();
    }
}