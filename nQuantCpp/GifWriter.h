#pragma once
#include "bitmapUtilities.h"

namespace GifEncode
{
	class GifWriter
	{
		private:
			bool _hasAlpha, _loop;
			wstring _destPath;
			unique_ptr<long[]> _frameDelays;
			unique_ptr<PropertyItem> _frameDelay, _loopPropertyItem;
		
			public:
				GifWriter(const wstring& destPath, const bool hasAlpha, const int count = 1, const long delay = 850, const bool loop = true);
				Status AddImages(vector<shared_ptr<Bitmap> >& bitmaps);
	};
}
