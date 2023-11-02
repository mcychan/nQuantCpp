/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2023 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

#include "stdafx.h"
#include "GilbertCurve.h"
#include "BlueNoise.h"
#include "CIELABConvertor.h"

#include <memory>
#include <list>
#include <algorithm>

namespace Peano
{
	struct ErrorBox
	{
		double yDiff = 0;
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

	bool sortedByYDiff;
	UINT m_width, m_height;
	const ARGB* m_image;
	const ColorPalette* m_pPalette;
	unsigned short* m_qPixels;
	DitherFn m_ditherFn;
	float* m_saliencies;
	GetColorIndexFn m_getColorIndexFn;
	list<ErrorBox> errorq;
	vector<float> m_weights;
	short* m_lookup;
	static BYTE DITHER_MAX = 9, ditherMax;
	static int margin, thresold;
	static const float BLOCK_SIZE = 343.0f;

	template <typename T> int sign(T val) {
		return (T(0) < val) - (val < T(0));
	}

	template<typename T>
	void insert_in_order(list<T>& items, T& element, int (*compareFn)(const T& a, const T& b)) {
		auto begin = items.begin();
		auto end = items.end();
		while (begin != end && compareFn(*begin, element) < 0)
			++begin;

		items.insert(begin, element);
	}

	int log2(int v) {
		int p = 0;
		int v2 = v >> 1;
		while (v2 > 0) {
			v2 >>= 1;
			++p;
		}		
		return p;
	}

	void initWeights(int size) {
		/* Dithers all pixels of the image in sequence using
		 * the Gilbert path, and distributes the error in
		 * a sequence of pixels size.
		 */
		errorq.resize(size);
		const auto weightRatio = (float)pow(BLOCK_SIZE + 1.0f, 1.0f / (size - 1.0f));
		auto weight = 1.0f;
		auto sumweight = 0.0f;
		m_weights.resize(size);
		for (int c = 0; c < size; ++c) {
			sumweight += (m_weights[size - c - 1] = 1.0f / weight);
			weight *= weightRatio;
		}

		weight = 0.0f; /* Normalize */
		for (int c = 0; c < size; ++c)
			weight += (m_weights[c] /= sumweight);
		m_weights[0] += 1.0f - weight;
	}

	inline int compare(const ErrorBox& o1, const ErrorBox& o2)
	{
		if (o1.yDiff > o2.yDiff)
			return -1;
		return o1.yDiff < o2.yDiff ? 1 : 0;
	}

	void ditherPixel(int x, int y)
	{
		int bidx = x + y * m_width;
		Color pixel(m_image[bidx]);
		ErrorBox error(pixel);
		int i = sortedByYDiff ? m_weights.size() - 1 : 0;
		auto maxErr = DITHER_MAX - 1;
		for (auto& eb : errorq) {
			if (i < 0 || i >= m_weights.size())
				break;

			for (int j = 0; j < eb.length(); ++j) {
				error[j] += eb[j] * m_weights[i];
				if (error[j] > maxErr)
					maxErr = error[j];
			}
			i += sortedByYDiff ? -1 : 1;
		}

		auto r_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[0], 0)));
		auto g_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[1], 0)));
		auto b_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[2], 0)));
		auto a_pix = static_cast<BYTE>(min(BYTE_MAX, max(error[3], 0)));

		Color c2 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
		if (m_pPalette->Count <= 32 && a_pix > 0xF0)
		{
			int offset = m_getColorIndexFn(c2);
			if (!m_lookup[offset])
				m_lookup[offset] = m_ditherFn(m_pPalette, c2.GetValue(), bidx) + 1;
			m_qPixels[bidx] = m_lookup[offset] - 1;

			if (m_saliencies != nullptr && CIELABConvertor::Y_Diff(pixel, c2) > margin) {
				auto strength = 1 / 3.0f;
				c2 = BlueNoise::diffuse(pixel, m_pPalette->Entries[m_qPixels[bidx]], 1.0f / m_saliencies[bidx], strength, x, y);
				m_qPixels[bidx] = m_ditherFn(m_pPalette, c2.GetValue(), bidx);
			}
		}
		else
			m_qPixels[bidx] = m_ditherFn(m_pPalette, c2.GetValue(), bidx);

		if (errorq.size() >= DITHER_MAX)
			errorq.pop_front();
		else if (!errorq.empty())
			initWeights(errorq.size());

		c2 = m_pPalette->Entries[m_qPixels[bidx]];
		error[0] = r_pix - c2.GetR();
		error[1] = g_pix - c2.GetG();
		error[2] = b_pix - c2.GetB();
		error[3] = a_pix - c2.GetA();

		auto denoise = m_pPalette->Count > 2;
		auto diffuse = BlueNoise::TELL_BLUE_NOISE[bidx & 4095] > thresold;
		auto illusion = !diffuse && BlueNoise::TELL_BLUE_NOISE[(int)(error.yDiff * 4096) & 4095] > thresold;
		error.yDiff = sortedByYDiff ? CIELABConvertor::Y_Diff(pixel, c2) : 1;

		int errLength = denoise ? error.length() - 1 : 0;
		for (int j = 0; j < errLength; ++j) {
			if (abs(error.p[j]) >= ditherMax) {
				if (diffuse)
					error[j] = (float)tanh(error.p[j] / maxErr * 8) * (ditherMax - 1);
				else if (illusion)
					error[j] = (float)(error.p[j] / maxErr * error.yDiff) * (ditherMax - 1);
				else
					error[j] /= (float)(1 + _sqrt(ditherMax));
			}
		}

		if (sortedByYDiff)
			insert_in_order<ErrorBox>(errorq, error, &compare);
		else
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
			for (int i = 0; i < w; ++i) {
				ditherPixel(x, y);
				x += dax;
				y += day;
			}
			return;
		}

		if (w == 1) {
			for (int i = 0; i < h; ++i) {
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

	void GilbertCurve::dither(const UINT width, const UINT height, const ARGB* pixels, const ColorPalette* pPalette, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, float* saliencies, double weight)
	{
		m_width = width;
		m_height = height;
		m_image = pixels;
		m_pPalette = pPalette;
		m_qPixels = qPixels;
		m_ditherFn = ditherFn;
		m_saliencies = saliencies;
		m_getColorIndexFn = getColorIndexFn;
		auto hasAlpha = weight < 0;

		errorq.clear();
		weight = abs(weight);
		sortedByYDiff = !hasAlpha && pPalette->Count >= 128 && weight >= .04;
		DITHER_MAX = weight < .01 ? (weight > .0025) ? (BYTE)25 : 16 : 9;
		auto edge = hasAlpha ? 1 : exp(weight) + .25;
		ditherMax = (hasAlpha || DITHER_MAX > 9) ? (BYTE)sqr(_sqrt(DITHER_MAX) + edge) : DITHER_MAX;
		if (pPalette->Count / weight > 5000 && (weight > .045 || (weight > .01 && pPalette->Count <= 64)))
			ditherMax = (BYTE)sqr(5 + edge);
		else if (pPalette->Count / weight < 3200 && pPalette->Count > 16 && pPalette->Count < 256)
			ditherMax = (BYTE)sqr(5 + edge);
		margin = (int) sqr(log2(pPalette->Count) - 1);
		thresold = DITHER_MAX > 9 ? -112 : -64;
		auto pLookup = make_unique<short[]>(USHRT_MAX + 1);
		m_lookup = pLookup.get();

		if (!sortedByYDiff)
			initWeights(DITHER_MAX);

		if (width >= height)
			generate2d(0, 0, width, 0, 0, height);
		else
			generate2d(0, 0, 0, height, width, 0);
	}
}
