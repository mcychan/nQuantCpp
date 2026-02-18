/* Generalized Hilbert ("gilbert") space-filling curve for rectangular domains of arbitrary (non-power of two) sizes.
Copyright (c) 2021 - 2026 Miller Cy Chan
* A general rectangle with a known orientation is split into three regions ("up", "right", "down"), for which the function calls itself recursively, until a trivial path can be produced. */

#include "stdafx.h"
#include "GilbertCurve.h"
#include "BlueNoise.h"
#include "CIELABConvertor.h"

#include <list>

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

	bool m_hasAlpha, m_dither, sortedByYDiff;
	unsigned short m_nMaxColor;
	UINT m_width, m_height;
	float beta;
	double m_weight;
	const ARGB *m_image, *m_pPalette;
	unsigned short* m_qPixels;
	ARGB* m_qColorPixels;
	DitherFn m_ditherFn;
	float* m_saliencies;
	GetColorIndexFn m_getColorIndexFn;
	list<ErrorBox> errorq;
	vector<float> m_weights;
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
		return sign(o2.yDiff - o1.yDiff);
	}

	float normalDistribution(float x, float peak) {
		const float mean = .5f, stdDev = .1f;

		// Calculate the probability density function (PDF)
		auto exponent = -pow(x - mean, 2) / (2 * pow(stdDev, 2));
		auto pdf = (1 / (stdDev * _sqrt(2 * M_PI))) * exp(exponent);
		auto maxPdf = 1 / (stdDev * _sqrt(2 * M_PI)); // Peak at x = mean
		auto scaledPdf = (pdf / maxPdf) * peak;
		return (float) max(0.0, min(peak, scaledPdf));
	}

	unsigned short ditherPixel(int x, int y, Color c2, float beta)
	{
		int bidx = x + y * m_width;
		Color pixel(m_image[bidx]);
		auto r_pix = c2.GetR();
		auto g_pix = c2.GetG();
		auto b_pix = c2.GetB();
		auto a_pix = c2.GetA();

		auto qPixelIndex = m_qPixels[bidx];
		auto strength = 1 / 3.0f;
		int acceptedDiff = max(2, m_nMaxColor - margin);
		if (m_nMaxColor <= 4 && m_saliencies[bidx] > .2f && m_saliencies[bidx] < .25f)
			c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], beta * 2 / m_saliencies[bidx], strength, x, y);
		else if (m_nMaxColor <= 4 || CIELABConvertor::Y_Diff(pixel, c2) < (2 * acceptedDiff)) {
			if (m_nMaxColor <= 128 || BlueNoise::TELL_BLUE_NOISE[bidx & 4095] > 0) {
				if (m_nMaxColor > 64) {
					auto kappa = m_saliencies[bidx] < .6f ? beta * .15f / m_saliencies[bidx] : beta * .4f / m_saliencies[bidx];
					c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], kappa, strength, x, y);
				}
				else if (m_nMaxColor > 16 && m_weight < .005)
					c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], beta * normalDistribution(m_saliencies[bidx], .5f) + beta, strength, x, y);
				else
					c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], beta * .5f / m_saliencies[bidx], strength, x, y);
			}
		}

		auto gamma = (m_nMaxColor <= 32 && m_weight < .01 && m_weight > .007) ? 1 - beta : beta;
		if (m_nMaxColor > 4 && CIELABConvertor::Y_Diff(pixel, c2) > (gamma * acceptedDiff)) {
			if (margin > 6 || gamma > beta) {
				auto kappa = m_saliencies[bidx] < .4f ? beta * .4f * m_saliencies[bidx] : beta * .4f / m_saliencies[bidx];
				Color c1 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
				if (m_nMaxColor > 32 && m_saliencies[bidx] < .9)
					kappa = beta * normalDistribution(m_saliencies[bidx], 2.0f);
				else {
					if (m_weight >= .0015 && m_saliencies[bidx] < .6)
						c1 = pixel;
					if (m_weight < .005 && m_saliencies[bidx] < .6)
						kappa = beta * normalDistribution(m_saliencies[bidx], m_weight < .0008 ? 2.5f : 1.75f);
					else if (m_nMaxColor >= 32 || CIELABConvertor::Y_Diff(c1, c2) > (gamma * M_PI * acceptedDiff)) {
						auto ub = 1 - m_nMaxColor / 320.0;
						if (m_saliencies[bidx] > .15 && m_saliencies[bidx] < ub)
							kappa = beta * (!sortedByYDiff && m_weight < .0025 ? .55f : .5f) / m_saliencies[bidx];
						else
							kappa = beta * normalDistribution(m_saliencies[bidx], m_weight < .0025 ? 1.82f : 2.0f);
					}
				}

				c2 = BlueNoise::diffuse(c1, m_pPalette[qPixelIndex], kappa, strength, x, y);
			}
			else if (m_nMaxColor <= 32 && m_weight >= .004)
				c2 = BlueNoise::diffuse(c2, m_pPalette[qPixelIndex], beta * normalDistribution(m_saliencies[bidx], .25f), strength, x, y);
			else
				c2 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);
		}

		if (DITHER_MAX < 16 && m_nMaxColor > 4 && m_saliencies[bidx] < .6f && CIELABConvertor::Y_Diff(pixel, c2) > margin - 1)
			c2 = Color::MakeARGB(a_pix, r_pix, g_pix, b_pix);

		return m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);
	}

	void diffusePixel(int x, int y)
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
		auto qPixelIndex = m_qPixels[bidx];
		if (m_saliencies != nullptr && m_dither && !sortedByYDiff && (!m_hasAlpha || pixel.GetA() < a_pix)) {
			if (m_nMaxColor >= 256 && m_saliencies[bidx] > .99f)
				qPixelIndex = m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);
			else
				qPixelIndex = ditherPixel(x, y, c2, beta);
		}
		else if (m_nMaxColor <= 32 && a_pix > 0xF0)
		{
			qPixelIndex = m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);

			int acceptedDiff = max(2, m_nMaxColor - margin);
			if (m_saliencies != nullptr && (CIELABConvertor::Y_Diff(pixel, c2) > acceptedDiff || CIELABConvertor::U_Diff(pixel, c2) > (2 * acceptedDiff))) {
				auto strength = 1 / 3.0f;
				c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], 1 / m_saliencies[bidx], strength, x, y);
				qPixelIndex = m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);
			}
		}
		else
			qPixelIndex = m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);

		if (errorq.size() >= DITHER_MAX)
			errorq.pop_front();
		else if (!errorq.empty())
			initWeights(errorq.size());

		c2 = m_pPalette[qPixelIndex];
		if (m_qPixels)
			m_qPixels[bidx] = qPixelIndex;
		else if (m_hasAlpha)
			m_qColorPixels[bidx] = c2.GetValue();
		else {
			Color c0 = m_pPalette[0];
			m_qColorPixels[bidx] = GetARGBIndex(c2, false, c0.GetA() == 0);
		}

		error[0] = r_pix - c2.GetR();
		error[1] = g_pix - c2.GetG();
		error[2] = b_pix - c2.GetB();
		error[3] = a_pix - c2.GetA();

		auto denoise = m_nMaxColor > 2;
		auto diffuse = BlueNoise::TELL_BLUE_NOISE[bidx & 4095] > thresold;
		error.yDiff = sortedByYDiff ? CIELABConvertor::Y_Diff(pixel, c2) : 1;
		auto illusion = !diffuse && BlueNoise::TELL_BLUE_NOISE[(int)(error.yDiff * 4096) & 4095] > thresold;

		auto unaccepted = false;
		int errLength = denoise ? error.length() - 1 : 0;
		for (int j = 0; j < errLength; ++j) {
			if (abs(error.p[j]) >= ditherMax) {
				if (sortedByYDiff && m_saliencies != nullptr)
					unaccepted = true;

				if (diffuse)
					error[j] = (float)tanh(error.p[j] / maxErr * 20) * (ditherMax - 1);
				else if (illusion)
					error[j] = (float)(error.p[j] / maxErr * error.yDiff) * (ditherMax - 1);
				else
					error[j] /= (float)(1 + _sqrt(ditherMax));
			}

			if (sortedByYDiff && m_saliencies == nullptr && abs(error.p[j]) >= DITHER_MAX)
				unaccepted = true;
		}

		if (unaccepted) {
			if (m_saliencies != nullptr)
				qPixelIndex = ditherPixel(x, y, c2, beta);
			else if (CIELABConvertor::Y_Diff(pixel, c2) > 3 && CIELABConvertor::U_Diff(pixel, c2) > 3) {
				auto strength = 1 / 3.0f;
				c2 = BlueNoise::diffuse(pixel, m_pPalette[qPixelIndex], strength, strength, x, y);
				qPixelIndex = m_ditherFn(m_pPalette, m_nMaxColor, c2.GetValue(), bidx);
			}

			c2 = m_pPalette[qPixelIndex];
			if (m_qPixels)
				m_qPixels[bidx] = qPixelIndex;
			else if (m_hasAlpha)
				m_qColorPixels[bidx] = c2.GetValue();
			else {
				Color c0 = m_pPalette[0];
				m_qColorPixels[bidx] = GetARGBIndex(c2, false, c0.GetA() == 0);
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
				diffusePixel(x, y);
				x += dax;
				y += day;
			}
			return;
		}

		if (w == 1) {
			for (int i = 0; i < h; ++i) {
				diffusePixel(x, y);
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

	void doDither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColor, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, float* saliencies, double weight)
	{
		m_width = width;
		m_height = height;
		m_image = pixels;
		m_pPalette = pPalette;
		m_nMaxColor = nMaxColor;
		
		m_ditherFn = ditherFn;
		m_getColorIndexFn = getColorIndexFn;
		m_hasAlpha = weight < 0;
		m_saliencies = saliencies;

		errorq.clear();
		weight = m_weight = abs(weight);
		margin = weight < .0025 ? 12 : weight < .004 ? 8 : 6;
		sortedByYDiff = m_saliencies && m_nMaxColor >= 128 && weight >= .02 && (!m_hasAlpha || weight < .18);
		beta = m_nMaxColor > 4 ? (float) (.6f - .00625f * m_nMaxColor) : 1;
		if (m_nMaxColor > 4) {
			auto boundary = .005 - .0000625 * m_nMaxColor;
			beta = (float) (weight > boundary ? max(.25, beta - m_nMaxColor * weight) : min(1.5, beta + m_nMaxColor * weight));
			if (m_nMaxColor > 16 && m_nMaxColor <= 32 && weight < .003)
				beta += .075f;
			else if (weight < .0015 || (m_nMaxColor > 32 && m_nMaxColor < 256))
				beta += .1f;
			if (m_nMaxColor >= 64 && (weight > .012 && weight < .0125) || (weight > .025 && weight < .03))
				beta += .05f;
			else if (m_nMaxColor > 32 && m_nMaxColor < 64 && weight < .015)
				beta = .55f;
			else if (m_nMaxColor > 16 && m_nMaxColor <= 32 && weight >= .005)
				beta += .1f;
		}
		else
			beta *= .95f;

		if (m_nMaxColor > 64 || (m_nMaxColor > 4 && weight > .02))
			beta *= .4f;
		if (m_nMaxColor > 64 && weight < .02)
			beta = .18f;

		DITHER_MAX = weight < .015 ? (weight > .0025) ? (BYTE)25 : 16 : 9;
		if (weight > .99) {
			beta = weight;
			DITHER_MAX = 25;
		}

		auto edge = m_hasAlpha ? 1 : exp(weight) + .25;
		auto deviation = !m_hasAlpha && weight > .002 ? .25 : 1;
		ditherMax = (m_hasAlpha || DITHER_MAX > 9) ? (BYTE)sqr(_sqrt(DITHER_MAX) + edge * deviation) : (BYTE) (DITHER_MAX * 1.5);
		int density = m_nMaxColor > 16 ? 3200 : 1500;
		if (m_nMaxColor / weight > 5000 && (weight > .045 || (weight > .01 && m_nMaxColor < 64)))
			ditherMax = (BYTE)sqr(5 + edge);
		else if (weight < .03 && m_nMaxColor / weight < density && m_nMaxColor >= 16 && m_nMaxColor < 256)
			ditherMax = (BYTE)sqr(5 + edge);
		thresold = DITHER_MAX > 9 ? -112 : -64;

		if (!sortedByYDiff)
			initWeights(DITHER_MAX);

		if (width >= height)
			generate2d(0, 0, width, 0, 0, height);
		else
			generate2d(0, 0, 0, height, width, 0);
	}

	void GilbertCurve::dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColor, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, unsigned short* qPixels, float* saliencies, double weight, bool dither)
	{
		m_qPixels = qPixels;
		m_dither = dither;
		doDither(width, height, pixels, pPalette, nMaxColor, ditherFn, getColorIndexFn, saliencies, weight);
	}

	void GilbertCurve::dither(const UINT width, const UINT height, const ARGB* pixels, const ARGB* pPalette, const UINT nMaxColor, DitherFn ditherFn, GetColorIndexFn getColorIndexFn, ARGB* qPixels, float* saliencies, double weight, bool dither)
	{
		m_qColorPixels = qPixels;
		m_dither = dither;
		doDither(width, height, pixels, pPalette, nMaxColor, ditherFn, getColorIndexFn, saliencies, weight);
	}
}
