#pragma once
/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2022 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

#include "stdafx.h"
#include "GilbertCurve.h"

#include <memory>
#include <deque>

namespace Peano
{	
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
        inline BYTE length() const
        {
            return 4;
        }
    };
	
    float m_divisor = 1.0f;
    UINT m_width, m_height;
    const ARGB* m_image;
    const ColorPalette* m_pPalette;
    unsigned short* m_qPixels;
    DitherFn m_ditherFn;
    GetColorIndexFn m_getColorIndexFn;
    deque<ErrorBox> errorq;
    float* m_weights;
    
    static const BYTE DITHER_MAX = 9;
    static const float BLOCK_SIZE = 343.0f; 
    
    template <typename T> int sign(T val) {
        return (T(0) < val) - (val < T(0));
    }
	
    void ditherPixel(int x, int y)
    {
        int bidx = x + y * m_width;
        Color pixel(m_image[bidx]);
        ErrorBox error(pixel);
        int i = 0;
        for (auto& eb : errorq) {
		    for(int j = 0; j < eb.length(); ++j)
                error[j] += eb[j] * m_weights[i];
		    ++i;
        }

        auto r_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[0], 0)));
        auto g_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[1], 0)));
        auto b_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[2], 0)));
        auto a_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[3], 0)));
		
        Color c2 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
        m_qPixels[bidx] = m_ditherFn(m_pPalette, c2.GetValue(), bidx);

        errorq.pop_front();
        c2 = m_pPalette->Entries[m_qPixels[bidx]];
        error[0] = r_pix - c2.GetR();
        error[1] = g_pix - c2.GetG();
        error[2] = b_pix - c2.GetB();
        error[3] = a_pix - c2.GetA();
		
        if (m_divisor < 3 || m_pPalette->Count > 16) {
            for (int j = 0; j < error.length(); ++j) {
                if (abs(error[j]) < DITHER_MAX)
                    continue;

                error[j] /= m_divisor;
            }
        }
        errorq.emplace_back(error);
    }
    
    void generate2d(int x, int y, int ax, int ay, int bx, int by) {    	
    	int w = abs(ax + ay);
    	int h = abs(bx + by);
    	int dax = sign(ax);
    	int day = sign(ay);
    	int dbx = sign(bx);
    	int dby = sign(by);

    	if (h == 1) {
    		for (int i = 0; i < w; ++i){
    			ditherPixel(x, y);
    			x += dax;
    			y += day;
    		}
    		return;
    	}

    	if (w == 1) {
    		for (int i = 0; i < h; ++i){
    			ditherPixel(x, y);
    			x += dbx;
    			y += dby;
    		}
    		return;
    	}

    	int ax2 = ax / 2;
    	int ay2 = ay / 2;
    	int bx2 = bx / 2;
    	int by2 = by / 2;

    	int w2 = abs(ax2 + ay2);
    	int h2 = abs(bx2 + by2);

    	if (2 * w > 3 * h) {
    		if ((w2 % 2) != 0 && w > 2) {
    			ax2 += dax;
    			ay2 += day;
    		}    		
    		generate2d(x, y, ax2, ay2, bx, by);
    		generate2d(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
    		return;
    	}
    	
    	if ((h2 % 2) != 0 && h > 2) {
    		bx2 += dbx;
    		by2 += dby;
    	}
		
    	generate2d(x, y, bx2, by2, ax2, ay2);
    	generate2d(x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
    	generate2d(x + (ax - dax) + (bx2 - dbx), y + (ay - day) + (by2 - dby), -bx2, -by2, -(ax - ax2), -(ay - ay2));    		
    }
	
    void GilbertCurve::dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, float divisor)
    {
        m_divisor = (divisor < 3) ? 0.4f + divisor - pPalette->Count / 64.0f : divisor;
    	if (divisor >= 1.5f || m_divisor > 1.5f)
            m_divisor = divisor;
        m_width = width;
        m_height = height;
        m_image = pixels;
        m_pPalette = pPalette;
        m_qPixels = qPixels;
        m_ditherFn = ditherFn;
        m_getColorIndexFn = getColorIndexFn;
        auto pWeights = make_unique<float[]>(DITHER_MAX);
        m_weights = pWeights.get();
    	
        /* Dithers all pixels of the image in sequence using
         * the Gilbert path, and distributes the error in
         * a sequence of 9 pixels.
         */
        errorq.clear();
        errorq.resize(DITHER_MAX);
        const auto weightRatio = (float)pow(BLOCK_SIZE + 1.0f, 1.0f / (DITHER_MAX - 1.0f));
        auto weight = 1.0f, sumweight = 0.0f;
        for (int c = 0; c < DITHER_MAX; ++c) {            
            sumweight += (m_weights[DITHER_MAX - c - 1] = 1.0f / weight);
            weight *= weightRatio;
        }

        weight = 0.0f; /* Normalize */
        for (int c = 0; c < DITHER_MAX; ++c)
            weight += (m_weights[c] /= sumweight);
        m_weights[0] += 1.0f - weight;
        
    	if (width >= height)
    		generate2d(0, 0, width, 0, 0, height);
    	else
    		generate2d(0, 0, 0, height, width, 0);
    }
}
