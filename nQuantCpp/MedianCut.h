/*
Copyright(c) 1989-1991 Jef Poskanzer.
Copyright(c) 1997-2002 Greg Roelofs; based on an idea by Stefan Schneider.
Copyright(c) 2009-2015 by Kornel Lesiński.
Copyright(c) 2015 Hao-Zhi Huang
Copyright (c) 2018 Miller Cy Chan

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once
#include <string>
#include <limits>
#include "EdgeAwareSQuantizer.h"

using namespace EdgeAwareSQuant;

//#define VITER_CACHE_LINE_GAP ((64+sizeof(viter_state)-1)/sizeof(viter_state))

namespace MedianCutQuant
{
	class MedianCut
	{
	public:
		virtual int quantizeImg(const vector<ARGB>& pixels, const UINT& width, Mat<float>& saliencyMap_float, ColorPalette* pPalette, UINT& newcolors);
		bool QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither = true);
	};
}