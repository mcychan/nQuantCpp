#pragma once
#include "bitmapUtilities.h"

namespace GifEncode
{
	class GifWriter
	{
		private:
			bool _hasAlpha, _loop;
			wstring _destPath;
			long _delay;
		
			public:
				GifWriter(const wstring& destPath, const bool hasAlpha, const long delay = 850, const bool loop = true);
				Status AddImages(vector<shared_ptr<Bitmap> >& bitmaps);
	};
}
