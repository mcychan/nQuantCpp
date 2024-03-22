#include "stdafx.h"
#include "GifWriter.h"

namespace GifEncode
{
	const int UintBytes = sizeof(long);

	GifWriter::GifWriter(const wstring& destPath, const bool hasAlpha, const long delay, const bool loop)
	{
		_destPath = destPath;
		_hasAlpha = hasAlpha;
		_delay = delay;
		_loop = loop;
	}

	static GUID GetEncoder()
	{
		const CLSID gifEncoderClsId = { 0x557cf402, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
		return gifEncoderClsId;
	}

	static void SetFrameDelay(Bitmap* pFirstBitmap, PropertyItem& frameDelay, vector<long>& frameDelays)
	{
		// PropertyItem for the frame delay (apparently, no other way to create a fresh instance).
		frameDelay.id = PropertyTagFrameDelay;
		frameDelay.type = PropertyTagTypeLong;
		// Length of the value in bytes.
		frameDelay.length = frameDelays.size() * UintBytes;
		// The value is an array of 4-byte entries: one per frame.
		// Every entry is the frame delay in 1/100-s of a second, in little endian.
		// E.g., here, we're setting the delay of every frame to 1 second.
		frameDelay.value = frameDelays.data();
		pFirstBitmap->SetPropertyItem(&frameDelay);
	}

	static void SetLoop(Bitmap* pFirstBitmap, PropertyItem& loopPropertyItem, const bool loop)
	{
		if (!loop)
			return;

		loopPropertyItem.id = PropertyTagLoopCount;
		loopPropertyItem.type = PropertyTagTypeShort;
		loopPropertyItem.length = 2;
		// 0 means to animate forever.
		short sValue = 0;
		loopPropertyItem.value = &sValue;
		pFirstBitmap->SetPropertyItem(&loopPropertyItem);
	}

	Status GifWriter::AddImages(vector<shared_ptr<Bitmap> >& bitmaps)
	{
		auto& pFirstBitmap = bitmaps[0];
		auto gifGUID = GetEncoder();

		// Params of the first frame.
		EncoderParameters encoderParams;
		encoderParams.Count = 1;
		encoderParams.Parameter[0].Guid = EncoderSaveFlag;
		encoderParams.Parameter[0].NumberOfValues = 1;
		encoderParams.Parameter[0].Type = EncoderParameterValueTypeLong;
		auto valueMultiFrame = EncoderValueMultiFrame;
		encoderParams.Parameter[0].Value = &valueMultiFrame;

		PropertyItem frameDelay;
		vector<long> frameDelays(bitmaps.size(), _delay / 10);
		SetFrameDelay(pFirstBitmap.get(), frameDelay, frameDelays);
		PropertyItem loopPropertyItem;
		SetLoop(pFirstBitmap.get(), loopPropertyItem, _loop);
		pFirstBitmap->Save(_destPath.c_str(), &gifGUID, &encoderParams);

		// Params of other frames.
		auto valueFrameDimensionTime = EncoderValueFrameDimensionTime;
		encoderParams.Parameter[0].Value = &valueFrameDimensionTime;
		Status status;
		for (int i = 1; i < bitmaps.size(); ++i)
			status = pFirstBitmap->SaveAdd(bitmaps[i].get(), &encoderParams);

		// Params for the finalizing call.
		encoderParams.Parameter[0].Type = EncoderParameterValueTypeLong;
		auto valueFlush = EncoderValueFlush;
		encoderParams.Parameter[0].Value = &valueFlush;
		return pFirstBitmap->SaveAdd(&encoderParams);
	}
}
