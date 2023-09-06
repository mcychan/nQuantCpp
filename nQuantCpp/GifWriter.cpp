#include "stdafx.h"
#include "GifWriter.h"

#include <memory>

namespace GifEncode
{
	const int UintBytes = sizeof(int);

	GifWriter::GifWriter(const wstring& destPath, const bool hasAlpha, const int count, const long delay, const bool loop)
	{
		_destPath = destPath;
		_hasAlpha = hasAlpha;
		_loop = loop;

		// PropertyItem for the frame delay (apparently, no other way to create a fresh instance).
		_frameDelay = make_unique<PropertyItem>();
		_frameDelay->id = PropertyTagFrameDelay;
		_frameDelay->type = PropertyTagTypeLong;
		// Length of the value in bytes.
		_frameDelay->length = (count + 1) * UintBytes;

		// E.g., here, we're setting the delay of every frame.
		_frameDelays = make_unique<long[]>(count + 1);
		for (int j = 0; j < count + 1; ++j)
			_frameDelays[j] = delay / 10;
		_frameDelay->value = _frameDelays.get();

		// PropertyItem for the number of animation loops.
		_loopPropertyItem = make_unique<PropertyItem>();
		_loopPropertyItem->id = PropertyTagLoopCount;
		_loopPropertyItem->type = PropertyTagTypeShort;
		_loopPropertyItem->length = 2;
		// 0 means to animate forever.
		short sValue = 0;
		_loopPropertyItem->value = &sValue;
	}

	static GUID GetEncoder()
	{
		const CLSID gifEncoderClsId = { 0x557cf402, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
		return gifEncoderClsId;
	}

	Status GifWriter::AddImages(vector<shared_ptr<Bitmap> >& bitmaps)
	{
		auto firstBitmap = bitmaps[0];
		auto gifGUID = GetEncoder();

		// Params of the first frame.
		EncoderParameters encoderParams;
		encoderParams.Count = 1;
		encoderParams.Parameter[0].Guid = EncoderSaveFlag;
		encoderParams.Parameter[0].NumberOfValues = 1;
		encoderParams.Parameter[0].Type = EncoderParameterValueTypeLong;
		long firstValue = EncoderValueMultiFrame;
		encoderParams.Parameter[0].Value = &firstValue;
		firstBitmap->SetPropertyItem(_frameDelay.get());
		if (_loop)
			firstBitmap->SetPropertyItem(_loopPropertyItem.get());
		firstBitmap->Save(_destPath.c_str(), &gifGUID, &encoderParams);

		// Params of other frames.
		firstValue = EncoderValueFrameDimensionTime;
		encoderParams.Parameter[0].Value = &firstValue;
		Status status;
		for (int i = 1; i < bitmaps.size(); ++i)
			status = firstBitmap->SaveAdd(bitmaps[i].get(), &encoderParams);

		// Params for the finalizing call.
		encoderParams.Parameter[0].Type = EncoderParameterValueTypeLong;
		firstValue = EncoderValueFlush;
		encoderParams.Parameter[0].Value = &firstValue;
		return firstBitmap->SaveAdd(&encoderParams);
	}
}
