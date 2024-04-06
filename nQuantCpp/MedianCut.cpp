/*
Copyright(c) 1989-1991 Jef Poskanzer.
Copyright(c) 1997-2002 Greg Roelofs; based on an idea by Stefan Schneider.
Copyright(c) 2009-2015 by Kornel Lesiński.
Copyright(c) 2015 Hao-Zhi Huang
Copyright(c) 2018-2021 Miller Cy Chan

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

#include "stdafx.h"
#include "MedianCut.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include "BlueNoise.h"
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace MedianCutQuant
{
	double PR = .299, PG = .587, PB = .114;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;
	unordered_map<ARGB, vector<unsigned short> > closestMap;
	unordered_map<ARGB, unsigned short> nearestMap;

	void getLab(const Color& c, CIELABConvertor::Lab& lab1)
	{
		auto got = pixelMap.find(c.GetValue());
		if (got == pixelMap.end()) {
			CIELABConvertor::RGB2LAB(c, lab1);
			pixelMap[c.GetValue()] = lab1;
		}
		else
			lab1 = got->second;
	}

	struct FloatPixel {
		FloatPixel()
		{
			a = 0, r = 0, g = 0, b = 0;
		}
		FloatPixel(ARGB argb)
		{
			Color c(argb);
			a = c.GetA() / 255.f;
			r = c.GetR() / 255.f;
			g = c.GetG() / 255.f;
			b = c.GetB() / 255.f;

		}
		FloatPixel(float _a, float _r, float _g, float _b)
		{
			a = _a, r = _r, g = _g, b = _b;
		}
		float a = 0, r = 0, g = 0, b = 0;
	};

	struct ChannelVariance {
		UINT chan;
		float variance;
	};

	struct SortTmp {
		float radius;
		UINT index;
	};

	struct ViterState {
		double a = 0, r = 0, g = 0, b = 0, total = 0;
	};

	struct Head {
		// colors less than radius away from vantage_point color will have best match in candidates
		FloatPixel vantage_point = FloatPixel(0, 0, 0, 0);
		float radius = 0;
		UINT num_candidates = 0;
		unique_ptr<FloatPixel[]> candidates_color;
		unique_ptr<UINT[]> candidates_index;
	};

	class ColorMapItem {
	public:
		FloatPixel fcolor;
		double popularity = 0;
	};

	class ColorMap {
	public:
		ColorMap(UINT _colors) {
			reset(_colors);
		}

		void reset(UINT _colors) {
			colors = _colors;
			palette = make_unique<ColorMapItem[]>(colors);
		}
		UINT colors;
		shared_ptr<ColorMap> subset_palette;
		unique_ptr<ColorMapItem[]> palette;
	};

	struct NearestMap {
		const ColorMap* map;
		vector<float> nearest_other_color_dist;
		vector<Head> heads;
		void reset(ULONG heads_size) {
			heads.resize(heads_size);
		}
	};

	static const double internal_gamma = 0.5499;

	static void toFloatSetGamma(float* gamma_lut, const double gamma, const UINT nMaxColors)
	{
		for (UINT i = 0; i < nMaxColors; ++i)
			gamma_lut[i] = pow((double)i / 255.0, internal_gamma / gamma);
	}

	inline static void rgbaToFloat(float* gamma_lut, const ARGB pixel, FloatPixel& f)
	{
		const float colorImportance = 1.0f;
		Color c(pixel);
		f.r = gamma_lut[c.GetR()] * colorImportance;
		f.g = gamma_lut[c.GetG()] * colorImportance;
		f.b = gamma_lut[c.GetB()] * colorImportance;
		f.a = c.GetA() / 255.0 * colorImportance;
	}

	inline static void rgbaToFloat(float* gamma_lut, vector_fixed<float, 4> & pixel, FloatPixel& f)
	{
		const float colorImportance = 1.0f;
		f.r = gamma_lut[static_cast<int>(pixel[0])] * colorImportance;
		f.g = gamma_lut[static_cast<int>(pixel[1])] * colorImportance;
		f.b = gamma_lut[static_cast<int>(pixel[2])] * colorImportance;
		f.a = pixel[3] / 255.0 * colorImportance;
	}

	static void floatToRgba(float gamma, const FloatPixel& px, vector_fixed<float, 4>& pixel)
	{
		if (px.a == 0) {
			pixel[0] = pixel[1] = pixel[2] = pixel[3] = 0;
			return;
		}

		float r = px.r / px.a,
			g = px.g / px.a,
			b = px.b / px.a,
			a = px.a;

		if (gamma != internal_gamma) {
			r = powf(r, gamma / internal_gamma);
			g = powf(g, gamma / internal_gamma);
			b = powf(b, gamma / internal_gamma);
		}

		// 256, because numbers are in range 1..255.9999m rounded down
		r *= 256.f;
		g *= 256.f;
		b *= 256.f;
		a *= 256.f;

		pixel[0] = r >= 255.f ? 1.f : r;
		pixel[1] = g >= 255.f ? 1.f : g;
		pixel[2] = b >= 255.f ? 1.f : b;
		pixel[3] = a >= 255.f ? 1.f : a;
	}

	class Box {
	public:
		FloatPixel color;
		FloatPixel variance;
		double sum = 0, total_error = 0, max_error = 0;
		unsigned short ind = 0, colors = 0;
	};

	class ColorHistArrItem {
	public:
		ARGB color;
		float perceptual_weight = 0;
		vector<pair<int, int> > pixLocation;
	};

	class ColorHistArrHead {
	public:
		UINT used = 0, capacity = 0;
		ColorHistArrItem inline1, inline2;
		unique_ptr<ColorHistArrItem[]> other_items;
	};

	class ColorHashTable {
	public:
		UINT ignorebits, maxcolors, colors;
		UINT hash_size;
		unique_ptr<ColorHistArrHead[]> buckets;

		ColorHashTable(UINT _hash_size, UINT _maxcolors, UINT _ignorebits) {
			hash_size = _hash_size;
			maxcolors = _maxcolors;
			ignorebits = _ignorebits;
			colors = 0;
			buckets = make_unique<ColorHistArrHead[]>(hash_size);
		}

		bool computeColorHash(const vector<ARGB>& pixels, const UINT& width, Mat<float>& importanceMap) {
			const UINT height = pixels.size() / width;
			const UINT channel_mask = 255U >> ignorebits << ignorebits;
			const UINT channel_hmask = (255U >> ignorebits) ^ 0xFFU;
			const UINT posterize_mask = channel_mask << 24 | channel_mask << 16 | channel_mask << 8 | channel_mask;
			const UINT posterize_high_mask = channel_hmask << 24 | channel_hmask << 16 | channel_hmask << 8 | channel_hmask;

			/* Go through the entire image, building a hash table of colors. */
			UINT pixelIndex = 0;
			for (int row = 0; row < height; ++row) {
				for (int col = 0; col < width; ++col) {
					float boost = 0.5 + importanceMap(row, col);
					/*if (importance_map) {
						boost = 0.5f + (double)*importance_map++ / 255.f;
					}*/

					// RGBA color is casted to long for easier hasing/comparisons
					auto& color = pixels[pixelIndex++];
					UINT hash = color % hash_size;

					/* head of the hash function stores first 2 colors inline (achl->used = 1..2),
					to reduce number of allocations of achl->other_items.
					*/
					auto& hashHead = buckets[hash];
					if (hashHead.inline1.color == color && hashHead.used) {
						hashHead.inline1.perceptual_weight += boost;
						hashHead.inline1.pixLocation.emplace_back(row, col);
						continue;
					}
					if (hashHead.used) {
						if (hashHead.used > 1) {
							if (hashHead.inline2.color == color) {
								hashHead.inline2.perceptual_weight += boost;
								hashHead.inline2.pixLocation.emplace_back(row, col);
								continue;
							}
							// other items are stored as an array (which gets reallocated if needed)
							auto other_items = hashHead.other_items.get();
							UINT i = 0;
							for (; i < hashHead.used - 2; ++i) {
								if (other_items[i].color == color) {
									other_items[i].perceptual_weight += boost;
									other_items[i].pixLocation.emplace_back(row, col);
									goto continue_outer_loop;
								}
							}

							// the array was allocated with spare items
							if (i < hashHead.capacity) {
								other_items[i].color = color;
								other_items[i].perceptual_weight = boost;
								other_items[i].pixLocation.emplace_back(row, col);
								hashHead.used++;
								++colors;
								continue;
							}

							if (++colors > maxcolors)
								return false;

							unique_ptr<ColorHistArrItem[]> pNew_items;
							UINT capacity;
							if (!other_items) { // there was no array previously, alloc "small" array
								capacity = 8;
								pNew_items = make_unique<ColorHistArrItem[]>(capacity);
							}
							else {
								// simply reallocs and copies array to larger capacity
								capacity = hashHead.capacity * 2 + 16;
								pNew_items = make_unique<ColorHistArrItem[]>(capacity);
								for (int i = 0; i < hashHead.capacity; ++i)
									pNew_items[i] = other_items[i];
							}
							pNew_items[i].color = color;
							pNew_items[i].perceptual_weight = boost;
							pNew_items[i].pixLocation.emplace_back(row, col);

							hashHead.other_items = move(pNew_items);
							hashHead.capacity = capacity;
							hashHead.used++;
						}
						else {
							// these are elses for first checks whether first and second inline-stored colors are used
							hashHead.inline2.color = color;
							hashHead.inline2.perceptual_weight = boost;
							hashHead.inline2.pixLocation.emplace_back(row, col);
							hashHead.used = 2;
							++colors;
						}
					}
					else {
						hashHead.inline1.color = color;
						hashHead.inline1.perceptual_weight = boost;
						hashHead.inline1.pixLocation.emplace_back(row, col);
						hashHead.used = 1;
						++colors;
					}

				continue_outer_loop:;
				}

			}
			return true;
		}
	};

	class HistItem {
	public:
		FloatPixel fcolor;
		float adjusted_weight,   // perceptual weight changed to tweak how mediancut selects colors
			perceptual_weight; // number of pixels weighted by importance of different areas of the picture

		float color_weight;      // these two change every time histogram subset is sorted
		union {
			UINT sort_value;
			BYTE likely_colormap_index;
		} tmp;

		vector<pair<int, int> > pixLocation;

		HistItem()
		{
			fcolor.a = 0;
			fcolor.r = 0;
			fcolor.g = 0;
			fcolor.b = 0;
			adjusted_weight = 0;
			perceptual_weight = 0;
			color_weight = 0;
			tmp.sort_value = 0;
		}

		HistItem(const HistItem& ref)
		{
			fcolor = ref.fcolor;
			adjusted_weight = ref.adjusted_weight;
			perceptual_weight = ref.perceptual_weight;
			color_weight = ref.color_weight;
			tmp = ref.tmp;
			pixLocation = ref.pixLocation;
		}

		HistItem operator=(const HistItem& ref)
		{
			fcolor = ref.fcolor;
			adjusted_weight = ref.adjusted_weight;
			perceptual_weight = ref.perceptual_weight;
			color_weight = ref.color_weight;
			tmp = ref.tmp;
			pixLocation = ref.pixLocation;
			return *this;
		}
	};

	class Histogram {
	public:
		vector<HistItem> histIterms;
		float max_perceptual_weight;
		double total_perceptual_weight;
		UINT size;
		UINT ignorebits;

		Histogram(const vector<ARGB>& pixels, const UINT& width, Mat<float>& importanceMap) {
			UINT ignorebits = 0;
			UINT maxcolors = 1966080;

			unique_ptr<ColorHashTable> ht;
			do {
				const UINT numPixel = pixels.size();
				UINT estimated_colors = min(maxcolors, numPixel / (ignorebits + (numPixel > 512 * 512 ? 5 : 4)));
				UINT hash_size = estimated_colors < 66000 ? 6673 : (estimated_colors < 200000 ? 12011 : 24019);
				ht = make_unique<ColorHashTable>(hash_size, maxcolors, ignorebits);

				// histogram uses noise contrast map for importance. Color accuracy in noisy areas is not very important.
				// noise map does not include edges to avoid ruining anti-aliasing

				bool added_ok = ht->computeColorHash(pixels, width, importanceMap);
				if (added_ok)
					break;
				if (!added_ok)
					++ignorebits;

			} while (!ht.get());

			hashToHist(*ht);
		}

	private:
		void hashToHist(const ColorHashTable& ht) {
			max_perceptual_weight = 0;
			ignorebits = ht.ignorebits;
			size = ht.colors;
			histIterms.reserve(size);

			/* Limit perceptual weight to 1/10th of the image surface area to prevent
			a single color from dominating all others. */
			float perceptualWeightLimit = 10;
			double total_weight = 0;

			auto gamma_lut = make_unique<float[]>(size);
			toFloatSetGamma(gamma_lut.get(), 1 / 2.2f, size);
			auto buckets = ht.buckets.get();
			for (UINT i = 0; i < ht.hash_size; ++i) {
				auto& achl = buckets[i];
				if (achl.used) {
					HistItem hi;
					rgbaToFloat(gamma_lut.get(), achl.inline1.color, hi.fcolor);
					hi.adjusted_weight = hi.perceptual_weight = min(achl.inline1.perceptual_weight, perceptualWeightLimit);
					hi.pixLocation = achl.inline1.pixLocation;
					if (hi.perceptual_weight > max_perceptual_weight)
						max_perceptual_weight = hi.perceptual_weight;
					total_weight += hi.adjusted_weight;
					histIterms.emplace_back(hi);

					if (achl.used > 1) {
						HistItem hi;
						rgbaToFloat(gamma_lut.get(), achl.inline2.color, hi.fcolor);
						hi.adjusted_weight = hi.perceptual_weight = min(achl.inline2.perceptual_weight, perceptualWeightLimit);
						hi.pixLocation = achl.inline2.pixLocation;
						if (hi.perceptual_weight > max_perceptual_weight)
							max_perceptual_weight = hi.perceptual_weight;
						total_weight += hi.adjusted_weight;
						histIterms.emplace_back(hi);
						auto other_items = achl.other_items.get();
						for (UINT k = 0; k < achl.used - 2; ++k) {
							HistItem hi;
							rgbaToFloat(gamma_lut.get(), other_items[k].color, hi.fcolor);
							hi.adjusted_weight = hi.perceptual_weight = min(other_items[k].perceptual_weight, perceptualWeightLimit);
							hi.pixLocation = other_items[k].pixLocation;
							if (hi.perceptual_weight > max_perceptual_weight)
								max_perceptual_weight = hi.perceptual_weight;
							total_weight += hi.adjusted_weight;
							histIterms.emplace_back(hi);
						}
					}
				}
			}

			total_perceptual_weight = total_weight;
		}
	};

	inline double variance_diff(double val, const double good_enough)
	{
		val *= val;
		if (val < sqr(good_enough))
			return val * 0.25;
		return val;
	}

	inline float channelDiff(const float x, const float y, const float alphas)
	{
		// maximum of channel blended on white, and blended on black
		// premultiplied alpha and backgrounds 0/1 shorten the formula
		const float black = x - y, white = black + alphas;
		return sqr(black) + sqr(white);
	}
	inline float colorDiff(const FloatPixel& px, const FloatPixel& py)
	{
		const float alphas = py.a - px.a;
		return channelDiff(px.r, py.r, alphas) + channelDiff(px.g, py.g, alphas) + channelDiff(px.b, py.b, alphas);
	}

	void avgColor(UINT clrs, const Histogram& hist, UINT beginIdx, const FloatPixel& center, FloatPixel& pixel)
	{
		double r = 0, g = 0, b = 0, a = 0, sum = 0;

		// first find final opacity in order to blend colors at that opacity
		for (UINT i = 0; i < clrs; ++i) {
			const auto& histIterm = hist.histIterms[beginIdx + i];
			const auto& px = histIterm.fcolor;
			sum += histIterm.adjusted_weight;
		}

		sum = 0;
		// reverse iteration for cache locality with previous loop
		for (int i = clrs - 1; i >= 0; --i) {
			double weight = 1.f;
			auto& histIterm = hist.histIterms[beginIdx + i];
			auto& px = histIterm.fcolor;

			/* give more weight to colors that are further away from average
			this is intended to prevent desaturation of images and fading of whites
			*/
			weight += sqr(center.r - px.r);
			weight += sqr(center.g - px.g);
			weight += sqr(center.b - px.b);

			weight *= histIterm.adjusted_weight;
			sum += weight;

			r += px.r * weight;
			g += px.g * weight;
			b += px.b * weight;
			a += px.a * weight;
		}

		if (!sum)
			sum = 1;

		a /= sum;
		r /= sum;
		g /= sum;
		b /= sum;

		pixel.a = a, pixel.r = r, pixel.g = g, pixel.b = b;
	}

	double boxError(const Box& box, const HistItem* hi)
	{
		const auto& avg = box.color;
		double total_error = 0;
		for (unsigned short i = 0; i < box.colors; ++i) {
			const int ind = box.ind + i;
			total_error += (double) colorDiff(avg, hi[ind].fcolor) * hi[ind].perceptual_weight;
		}

		return total_error;
	}

	bool totalBoxErrorBelowTarget(double target_mse, Box* bv, unsigned short boxes, const Histogram& hist)
	{
		target_mse *= hist.total_perceptual_weight;
		double total_error = 0;

		for (unsigned short i = 0; i < boxes; ++i) {
			// error is (re)calculated lazily
			if (bv[i].total_error >= 0)
				total_error += bv[i].total_error;

			if (total_error > target_mse)
				return false;
		}

		for (unsigned short i = 0; i < boxes; ++i) {
			if (bv[i].total_error < 0) {
				bv[i].total_error = boxError(bv[i], hist.histIterms.data());
				total_error += bv[i].total_error;
			}
			if (total_error > target_mse)
				return false;
		}

		return true;
	}

	FloatPixel boxVariance(Histogram& hist, const Box& box)
	{
		const auto& mean = box.color;
		double variancea = 0, variancer = 0, varianceg = 0, varianceb = 0;

		for (unsigned short i = 0; i < box.colors; ++i) {
			const auto& histIterm = hist.histIterms[box.ind + i];
			const auto& px = histIterm.fcolor;
			double weight = histIterm.adjusted_weight;
			variancea += variance_diff(mean.a - px.a, 2.0 / 256.0) * weight;
			variancer += variance_diff(mean.r - px.r, 1.0 / 256.0) * weight;
			varianceg += variance_diff(mean.g - px.g, 1.0 / 256.0) * weight;
			varianceb += variance_diff(mean.b - px.b, 1.0 / 256.0) * weight;
		}

		return FloatPixel(variancea * (4.0 / 16.0),
			variancer * (7.0 / 16.0),
			varianceg * (9.0 / 16.0),
			varianceb * (5.0 / 16.0));
	}

	double boxMaxError(Histogram& hist, const Box& box)
	{
		const auto& mean = box.color;
		double max_error = 0;

		for (unsigned short i = 0; i < box.colors; ++i) {
			const int ind = box.ind + i;
			if (ind >= hist.histIterms.size())
				return max_error;

			const double diff = colorDiff(mean, hist.histIterms[box.ind + i].fcolor);
			if (diff > max_error)
				max_error = diff;
		}
		return max_error;
	}

	void setColormapFromBoxes(const ColorMap& map, const Box* bv, unsigned short boxes, const Histogram& hist)
	{
		/*
		** Ok, we've got enough boxes.  Now choose a representative color for
		** each box.  There are a number of possible ways to make this choice.
		** One would be to choose the center of the box; this ignores any structure
		** within the boxes.  Another method would be to average all the colors in
		** the box - this is the method specified in Heckbert's paper.
		*/

		for (unsigned short bi = 0; bi < boxes; ++bi) {
			map.palette[bi].fcolor = bv[bi].color;

			/* store total color popularity (perceptual_weight is approximation of it) */
			map.palette[bi].popularity = 0;
			for (UINT i = bv[bi].ind; i < bv[bi].ind + bv[bi].colors; ++i) {
				const auto& histIterm = hist.histIterms[i];
				map.palette[bi].popularity += histIterm.perceptual_weight;
			}
		}
	}

	int bestSplittableBox(const Box* bv, unsigned short boxes, const double max_mse)
	{
		int bi = -1;
		double maxsum = 0;
		for (unsigned short i = 0; i < boxes; ++i) {
			if (bv[i].colors < 2)
				continue;

			// looks only at max variance, because it's only going to split by it
			const double cv = max(bv[i].variance.r, max(bv[i].variance.g, bv[i].variance.b));
			double thissum = bv[i].sum * max(bv[i].variance.a, cv);

			if (bv[i].max_error > max_mse)
				thissum *= bv[i].max_error / max_mse;

			if (thissum > maxsum) {
				maxsum = thissum;
				bi = i;
			}
		}
		return bi;
	}

	/* increase histogram popularity by difference from the final color (this is used as part of feedback loop) */
	void adjustHistogram(Histogram& hist, const ColorMap& map, const Box* bv, unsigned short boxes)
	{
		for (unsigned short bi = 0; bi < boxes; ++bi) {
			for (UINT i = bv[bi].ind; i < bv[bi].ind + bv[bi].colors; ++i) {
				auto& histIterm = hist.histIterms[i];
				auto& pixel = map.palette[bi];
				histIterm.adjusted_weight *= _sqrt(1.0 + colorDiff(pixel.fcolor, histIterm.fcolor) / 4.0);
				histIterm.tmp.likely_colormap_index = bi;
			}
		}
	}

	inline void histItemSwap(HistItem* l, HistItem* r)
	{
		if (l != r)
			swap(*l, *r);
	}

	inline UINT qsortPivot(const HistItem* const base, const UINT len)
	{
		if (len < 32)
			return len / 2;

		const UINT aidx = 8, bidx = len / 2, cidx = len - 1;
		const UINT a = base[aidx].tmp.sort_value, b = base[bidx].tmp.sort_value, c = base[cidx].tmp.sort_value;
		return (a < b) ? ((b < c) ? bidx : ((a < c) ? cidx : aidx))
			: ((b > c) ? bidx : ((a < c) ? aidx : cidx));
	}

	inline UINT qsortPartition(HistItem* const base, const UINT len)
	{
		UINT l = 1, r = len;
		if (len >= 8) {
			histItemSwap(&base[0], &base[qsortPivot(base, len)]);
		}

		const UINT pivot_value = base[0].tmp.sort_value;
		while (l < r) {
			if (base[l].tmp.sort_value >= pivot_value)
				l++;
			else {
				while (l < --r && base[r].tmp.sort_value <= pivot_value) {}
				histItemSwap(&base[l], &base[r]);
			}
		}
		histItemSwap(&base[0], &base[--l]);

		return l;
	}

	HistItem* histItemSortHalfVar(HistItem* base, UINT len, double* const lowervar, const double halfvar)
	{
		for (;;) {
			const UINT l = qsortPartition(base, len), r = l + 1;

			// check if sum of left side is smaller than half,
			// if it is, then it doesn't need to be sorted
			UINT t = 0;
			double tmpsum = *lowervar;
			while (t <= l && tmpsum < halfvar) {
				tmpsum += base[t].color_weight;
				++t;
			}

			if (tmpsum < halfvar)
				* lowervar = tmpsum;
			else {
				if (l > 0) {
					auto res = histItemSortHalfVar(base, l, lowervar, halfvar);
					if (res)
						return res;
				}
				else {
					// End of left recursion. This will be executed in order from the first element.
					*lowervar += base[0].color_weight;
					if (*lowervar > halfvar)
						return &base[0];
				}
			}

			if (len > r) {
				base += r;
				len -= r; // tail-recursive "call"
			}
			else {
				*lowervar += base[r].color_weight;
				return (*lowervar > halfvar) ? &base[r] : nullptr;
			}
		}
	}

	static int compareVariance(const void* ch1, const void* ch2)
	{
		return ((const ChannelVariance*)ch1)->variance > ((const ChannelVariance*)ch2)->variance ? -1 :
			(((const ChannelVariance*)ch1)->variance < ((const ChannelVariance*)ch2)->variance ? 1 : 0);
	}

	static int compareradius(const void* ap, const void* bp)
	{
		float a = ((const SortTmp*)ap)->radius;
		float b = ((const SortTmp*)bp)->radius;
		return a > b ? 1 : (a < b ? -1 : 0);
	}

	/** quick select algorithm */
	void histItemSortRange(HistItem* base, UINT len, UINT sort_start)
	{
		for (;;) {
			const UINT l = qsortPartition(base, len), r = l + 1;

			if (l > 0 && sort_start < l)
				len = l;
			else if (r < len && sort_start > r) {
				base += r;
				len -= r;
				sort_start -= r;
			}
			else
				break;
		}
	}

	/** finds median in unsorted set by sorting only minimum required */
	void getMedian(const Box& b, Histogram& hist, FloatPixel& median)
	{
		const UINT median_start = (b.colors - 1) / 2;

		histItemSortRange(&(hist.histIterms[b.ind]), b.colors, median_start);

		if (b.colors & 1)
			median = hist.histIterms[b.ind + median_start].fcolor;

		// technically the second color is not guaranteed to be sorted correctly
		// but most of the time it is good enough to be useful
		avgColor(2, hist, b.ind + median_start, FloatPixel(0.5, 0.5, 0.5, 0.5), median);
	}

	inline double colorWeight(const FloatPixel& median, const HistItem& h)
	{
		float diff = colorDiff(median, h.fcolor);
		// if color is "good enough", don't split further
		if (diff < 2.f / 256.f / 256.f)
			diff /= 2.f;
		return _sqrt(diff) * (_sqrt(1.0 + h.adjusted_weight) - 1.0);
	}

	inline float minChannelDiff(const float x, const float y, const float alphas)
	{
		const float black = x - y, white = black + alphas;
		return min(sqr(black), sqr(white)) * 2.f;
	}

	inline float minColorDiff(const FloatPixel& px, const FloatPixel& py)
	{
		const float alphas = py.a - px.a;
		return minChannelDiff(px.r, py.r, alphas) +
			minChannelDiff(px.g, py.g, alphas) +
			minChannelDiff(px.b, py.b, alphas);
	}

	double prepareSort(const Box& b, Histogram& hist)
	{
		/*
		** Sort dimensions by their variance, and then sort colors first by dimension with highest variance
		*/
		ChannelVariance channels[4] = {
			{ 1, b.variance.r },
		{ 2, b.variance.g },
		{ 3, b.variance.b },
		{ 0, b.variance.a },
		};

		qsort(channels, 4, sizeof(channels[0]), compareVariance);

		for (unsigned short i = 0; i < b.colors; ++i) {
			const float* chans = (const float*)& hist.histIterms[b.ind + i].fcolor;
			// Only the first channel really matters. When trying median cut many times
			// with different histogram weights, I don't want sort randomness to influence outcome.
			hist.histIterms[b.ind + i].tmp.sort_value = ((UINT)(chans[channels[0].chan] * 65535.0) << 16) |
				(UINT)((chans[channels[2].chan] + chans[channels[1].chan] / 2.0 + chans[channels[3].chan] / 4.0) * 65535.0);
		}

		FloatPixel median;
		getMedian(b, hist, median);

		// box will be split to make color_weight of each side even
		const unsigned short ind = b.ind, end = ind + b.colors;
		double totalvar = 0;
		for (unsigned short j = ind; j < end; ++j)
			totalvar += (hist.histIterms[j].color_weight = colorWeight(median, hist.histIterms[j]));
		return totalvar / 2.0;
	}

	unique_ptr<ColorMap> medianCut(Histogram& hist, const UINT& newcolors, const double target_mse, const double max_mse)
	{
		auto bv = make_unique<Box[]>(newcolors);

		/*
		** Set up the initial box.
		*/
		bv[0].ind = 0;
		bv[0].colors = hist.size;
		avgColor(bv[0].colors, hist, bv[0].ind, FloatPixel(0.5, 0.5, 0.5, 0.5), bv[0].color);
		bv[0].variance = boxVariance(hist, bv[0]);
		bv[0].max_error = boxMaxError(hist, bv[0]);
		bv[0].sum = 0;
		bv[0].total_error = -1;
		for (unsigned short i = 0; i < bv[0].colors; ++i)
			bv[0].sum += hist.histIterms[i].adjusted_weight;

		UINT boxes = 1;

		// remember smaller palette for fast searching
		unique_ptr<ColorMap> representative_subset;
		UINT subset_size = ceilf(powf(newcolors, 0.7f));

		/*
		** Main loop: split boxes until we have enough.
		*/
		while (boxes < newcolors) {

			if (boxes == subset_size) {
				representative_subset = make_unique<ColorMap>(boxes);
				setColormapFromBoxes(*representative_subset, bv.get(), boxes, hist);
			}

			// first splits boxes that exceed quality limit (to have colors for things like odd green pixel),
			// later raises the limit to allow large smooth areas/gradients get colors.
			const double current_max_mse = max_mse + (boxes / (double)newcolors) * _sqrt(newcolors) * max_mse;
			const int bi = bestSplittableBox(bv.get(), boxes, current_max_mse);
			if (bi < 0)
				break;        /* ran out of colors! */

			auto indx = bv[bi].ind;
			auto clrs = bv[bi].colors;

			/*
			Classic implementation tries to get even number of colors or pixels in each subdivision.

			Here, instead of popularity I use (sqrt(popularity)*variance) metric.
			Each subdivision balances number of pixels (popular colors) and low variance -
			boxes can be large if they have similar colors. Later boxes with high variance
			will be more likely to be split.

			Median used as expected value gives much better results than mean.
			*/

			const double halfvar = prepareSort(bv[bi], hist);
			double lowervar = 0;

			// hist_item_sort_halfvar sorts and sums lowervar at the same time
			// returns item to break at mminus one, which does smell like an off-by-one error.
			auto break_p = histItemSortHalfVar(&hist.histIterms[indx], clrs, &lowervar, halfvar);
			int break_at = clrs - 1;
			if (break_p)
				break_at = min(break_at, break_p - &hist.histIterms[indx] + 1);

			/*
			** Split the box.
			*/
			double sm = bv[bi].sum;
			double lowersum = 0;
			for (int i = 0; i < break_at; ++i)
				lowersum += hist.histIterms[indx + i].adjusted_weight;

			const auto& previous_center = bv[bi].color;
			bv[bi].colors = break_at;
			bv[bi].sum = lowersum;
			avgColor(bv[bi].colors, hist, bv[bi].ind, previous_center, bv[bi].color);
			bv[bi].total_error = -1;
			bv[bi].variance = boxVariance(hist, bv[bi]);
			bv[bi].max_error = boxMaxError(hist, bv[bi]);
			bv[boxes].ind = indx + break_at;
			bv[boxes].colors = clrs - break_at;
			bv[boxes].sum = sm - lowersum;
			avgColor(bv[boxes].colors, hist, bv[boxes].ind, previous_center, bv[boxes].color);
			bv[boxes].total_error = -1;
			bv[boxes].variance = boxVariance(hist, bv[boxes]);
			bv[boxes++].max_error = boxMaxError(hist, bv[boxes]);

			if (totalBoxErrorBelowTarget(target_mse, bv.get(), boxes, hist))
				break;
		}

		auto map = make_unique<ColorMap>(boxes);
		setColormapFromBoxes(*map, bv.get(), boxes, hist);
		map->subset_palette = move(representative_subset);
		adjustHistogram(hist, *map, bv.get(), boxes);
		return map;
	}

	void viterFinalize(ColorMap& map, ViterState* average_color)
	{
		for (UINT i = 0; i < map.colors; ++i) {
			double a = 0, r = 0, g = 0, b = 0, total = 0;

			// Aggregate results from all threads
			//for (UINT t = 0; t < max_threads; t++) {
			const UINT offset = i;

			a += average_color[offset].a;
			r += average_color[offset].r;
			g += average_color[offset].g;
			b += average_color[offset].b;
			total += average_color[offset].total;

			if (total)
				map.palette[i].fcolor = FloatPixel(a / total, r / total, g / total, b / total);

			map.palette[i].popularity = total;
		}
	}

	void viterUpdateColor(const FloatPixel& acolor, const float value, const ColorMap& map, UINT match, ViterState* average_color)
	{
		//match += thread * (VITER_CACHE_LINE_GAP + map.colors);
		const auto val = (double) value;
		average_color[match].a += acolor.a * val;
		average_color[match].r += acolor.r * val;
		average_color[match].g += acolor.g * val;
		average_color[match].b += acolor.b * val;
		average_color[match].total += val;
	}

	UINT findSlow(const FloatPixel& px, const ColorMap& map)
	{
		UINT best = 0;
		float bestdiff = colorDiff(px, map.palette[0].fcolor);

		for (UINT i = 1; i < map.colors; ++i) {
			float diff = colorDiff(px, map.palette[i].fcolor);
			if (diff < bestdiff) {
				best = i;
				bestdiff = diff;
			}
		}
		return best;
	}

	shared_ptr<ColorMap> getSubsetPalette(const ColorMap& map)
	{
		if (map.subset_palette.get())
			return map.subset_palette;

		UINT subset_size = (map.colors + 3) / 4;
		auto subset_palette = make_shared<ColorMap>(subset_size);

		for (UINT i = 0; i < subset_size; ++i)
			subset_palette->palette[i] = map.palette[i];
		return subset_palette;
	}

	void buildHead(Head& h, const FloatPixel& px, const ColorMap& map, UINT num_candidates, float error_margin, bool* skip_index, UINT* skipped)
	{
		auto colors = make_unique<SortTmp[]>(map.colors);
		UINT colorsused = 0;

		for (UINT i = 0; i < map.colors; ++i) {
			if (skip_index[i])
				continue; // colors in skip_index have been eliminated already in previous heads
			colors[colorsused].index = i;
			colors[colorsused++].radius = colorDiff(px, map.palette[i].fcolor);
		}

		qsort(colors.get(), colorsused, sizeof(colors[0]), compareradius);
		// ASSERT(colorsused < 2 || colors[0].radius <= colors[1].radius); // closest first

		num_candidates = min(colorsused, num_candidates);

		h.candidates_color = make_unique<FloatPixel[]>(num_candidates);
		h.candidates_index = make_unique<UINT[]>(num_candidates);
		h.vantage_point = FloatPixel(px.a, px.r, px.g, px.b);
		h.num_candidates = num_candidates;
		auto candidates_color = h.candidates_color.get();
		auto candidates_index = h.candidates_index.get();
		for (UINT i = 0; i < num_candidates; ++i) {
			candidates_color[i] = map.palette[colors[i].index].fcolor;
			candidates_index[i] = colors[i].index;
		}
		// if all colors within this radius are included in candidates, then there cannot be any other better match
		// farther away from the vantage point than half of the radius. Due to alpha channel must assume pessimistic radius.
		h.radius = minColorDiff(px, candidates_color[num_candidates - 1]) / 4.0f; // /4 = half of radius, but radius is squared

		for (UINT i = 0; i < num_candidates; ++i) {
			// divide again as that's matching certain subset within radius-limited subset
			// - 1/256 is a tolerance for miscalculation (seems like colordifference isn't exact)
			if (colors[i].radius < h.radius / 4.f - error_margin) {
				skip_index[colors[i].index] = true;
				(*skipped)++;
			}
		}
	}

	float distanceFromNearestOtherColor(const ColorMap& map, const UINT i)
	{
		float second_best = INT_MAX;
		for (UINT j = 0; j < map.colors; ++j) {
			if (i == j)
				continue;
			float diff = colorDiff(map.palette[i].fcolor, map.palette[j].fcolor);
			if (diff <= second_best)
				second_best = diff;
		}
		return second_best;
	}

	void nearestInit(const ColorMap& map, NearestMap& centroids)
	{
		auto subset_palette = getSubsetPalette(map);
		const UINT num_vantage_points = map.colors > 16 ? min(map.colors / 4, subset_palette->colors) : 0;
		const ULONG heads_size = num_vantage_points + 1; // +1 is fallback head

		centroids.nearest_other_color_dist.resize(map.colors);
		centroids.reset(heads_size);

		for (UINT i = 0; i < map.colors; ++i) {
			const float dist = distanceFromNearestOtherColor(map, i);
			centroids.nearest_other_color_dist[i] = dist / 4.f; // half of squared distance
		}

		centroids.map = &map;

		UINT skipped = 0;
		auto skip_index = make_unique<bool[]>(map.colors);
		for (int i = 0; i < map.colors; ++i)
			skip_index[i] = false;

		// floats and colordifference calculations are not perfect
		const float error_margin = 8.f / 256.f / 256.f;
		UINT h = 0;
		for (; h < num_vantage_points; ++h) {
			UINT num_candiadtes = 1 + (map.colors - skipped) / ((1 + num_vantage_points - h) / 2);

			buildHead(centroids.heads[h], subset_palette->palette[h].fcolor, map, num_candiadtes, error_margin, skip_index.get(), &skipped);
			if (centroids.heads[h].num_candidates == 0)
				break;
		}

		// assumption that there is no better color within radius of vantage point color
		// holds true only for colors within convex hull formed by palette colors.
		// since finding proper convex hull is more than a few lines, this
		// is a cheap shot at finding just few key points.
		const FloatPixel extrema[] = {
			FloatPixel(0, 0, 0, 0),

			FloatPixel(.5, 0, 0, 0), FloatPixel(.5, 1, 0, 0),
			FloatPixel(.5, 0, 0, 1), FloatPixel(.5, 1, 0, 1),
			FloatPixel(.5, 0, 1, 0), FloatPixel(.5, 1, 1, 0),
			FloatPixel(.5, 0, 1, 1), FloatPixel(.5, 1, 1, 1),

			FloatPixel(1, 0, 0, 0), FloatPixel(1, 1, 0, 0),
			FloatPixel(1, 0, 0, 1), FloatPixel(1, 1, 0, 1),
			FloatPixel(1, 0, 1, 0), FloatPixel(1, 1, 1, 0),
			FloatPixel(1, 0, 1, 1), FloatPixel(1, 1, 1, 1),

			FloatPixel(1, .5, 0, 0), FloatPixel(1, 0, .5, 0), FloatPixel(1, 0, 0, .5),
			FloatPixel(1, .5, 0, 1), FloatPixel(1, 0, .5, 1), FloatPixel(1, 0, 1, .5),
			FloatPixel(1, .5, 1, 0), FloatPixel(1, 1, .5, 0), FloatPixel(1, 1, 0, .5),
			FloatPixel(1, .5, 1, 1), FloatPixel(1, 1, .5, 1), FloatPixel(1, 1, 1, .5),
		};
		for (UINT i = 0; i < sizeof(extrema) / sizeof(extrema[0]); ++i)
			skip_index[findSlow(extrema[i], map)] = 0;

		buildHead(centroids.heads[h], FloatPixel(0, 0, 0, 0), map, map.colors, error_margin, skip_index.get(), &skipped);
		centroids.heads[h].radius = INT_MAX;
	}

	void adjustHistogramCallback(HistItem* item, float diff)
	{
		item->adjusted_weight = (item->perceptual_weight + item->adjusted_weight) * (sqrtf(1.f + diff));
	}

	UINT nearestSearch(NearestMap& centroids, const FloatPixel& px, int likely_colormap_index, float& diff)
	{
		// ASSERT(likely_colormap_index < centroids->map->colors);
		auto palette = centroids.map->palette.get();
		float guess_diff = colorDiff(palette[likely_colormap_index].fcolor, px);
		if (guess_diff < centroids.nearest_other_color_dist[likely_colormap_index]) {
			diff = guess_diff;
			return likely_colormap_index;
		}

		for (const auto& head : centroids.heads) {
			float vantage_point_dist = colorDiff(px, head.vantage_point);
			auto candidates_color = head.candidates_color.get();
			auto candidates_index = head.candidates_index.get();
			if (vantage_point_dist <= head.radius) {
				// ASSERT(heads[i].num_candidates);
				UINT ind = 0;
				float dist = colorDiff(px, candidates_color[0]);

				for (UINT j = 1; j < head.num_candidates; ++j) {
					float newdist = colorDiff(px, candidates_color[j]);

					if (newdist < dist) {
						dist = newdist;
						ind = j;
					}
				}
				if (diff)
					diff = dist;
				return candidates_index[ind];
			}
		}
		return likely_colormap_index;
	}

	double viterDoIteration(Histogram& hist, ColorMap& map, const bool first_run)
	{
		const UINT max_threads = omp_get_max_threads();
		auto average_color = make_unique<ViterState[]>(map.colors);
		NearestMap n;
		nearestInit(map, n);
		auto achv = hist.histIterms.data();
		const int hist_size = hist.size;

		double total_diff = 0;
		//#pragma omp parallel for if (hist_size > 3000) \
		//        schedule(static) default(none) shared(average_color,callback) reduction(+:total_diff)
		for (int j = 0; j < hist_size; ++j) {
			float diff = 0;
			UINT match = nearestSearch(n, achv[j].fcolor, achv[j].tmp.likely_colormap_index, diff);
			achv[j].tmp.likely_colormap_index = match; // which centroid it belongs to
			total_diff += (double) diff * achv[j].perceptual_weight;

			viterUpdateColor(achv[j].fcolor, achv[j].perceptual_weight, map, match, average_color.get());

			if (!first_run)
				adjustHistogramCallback(&achv[j], diff);
		}

		viterFinalize(map, average_color.get());
		return total_diff / hist.total_perceptual_weight;
	}

	int iterMedianCut(unique_ptr<ColorMap>& pFinalColorMap, Histogram& hist, UINT& newcolors, double& palette_error)
	{
		// if output is posterized it doesn't make sense to aim for perfrect colors, so increase target_mse
		// at this point actual gamma is not set, so very conservative posterization estimate is used
		const double target_mse = sqr(1 / 1024.0);
		double least_error = INT_MAX;
		const double max_mse = max(max(90.0 / 65536.0, target_mse), least_error) * 1.2;
		int feedback_loop_trials = 29;
		double target_mse_overshoot = feedback_loop_trials > 0 ? 1.05 : 1.0;
		const double percent = (double)(feedback_loop_trials > 0 ? feedback_loop_trials : 1) / 100.0;

		do {
			auto pNewMap = medianCut(hist, newcolors, target_mse * target_mse_overshoot, max_mse);

			// after palette has been created, total error (MSE) is calculated to keep the best palette
			// at the same time Voronoi iteration is done to improve the palette
			// and histogram weights are adjusted based on remapping error to give more weight to poorly matched colors
			const bool first_run_of_target_mse = !pFinalColorMap.get() && target_mse > 0;
			double total_error = viterDoIteration(hist, *pNewMap, first_run_of_target_mse);

			// goal is to increase quality or to reduce number of colors used if quality is good enough
			if (!pFinalColorMap.get() || total_error < least_error || (total_error <= target_mse && pNewMap->colors < newcolors)) {
				pFinalColorMap = move(pNewMap);

				if (total_error < target_mse && total_error > 0) {
					// voronoi iteration improves quality above what mediancut aims for
					// this compensates for it, making mediancut aim for worse
					target_mse_overshoot = min(target_mse_overshoot * 1.25, target_mse / total_error);
				}

				least_error = total_error;

				// if number of colors could be reduced, try to keep it that way
				// but allow extra color as a bit of wiggle room in case quality can be improved too
				newcolors = min(pFinalColorMap->colors + 1, newcolors);

				--feedback_loop_trials; // asymptotic improvement could make it go on forever
			}
			else {
				for (UINT j = 0; j < hist.size; ++j)
					hist.histIterms[j].adjusted_weight = (hist.histIterms[j].perceptual_weight + hist.histIterms[j].adjusted_weight) / 2.0;

				target_mse_overshoot = 1.0;
				feedback_loop_trials -= 6;
				// if error is really bad, it's unlikely to improve, so end sooner
				if (total_error > least_error * 4)
					feedback_loop_trials -= 3;
			}

		} while (feedback_loop_trials > 0);

		// likely_colormap_index (used and set in viter_do_iteration) can't point to index outside colormap
		if (pFinalColorMap->colors < 256) {
			for (UINT j = 0; j < hist.size; ++j) {
				if (hist.histIterms[j].tmp.likely_colormap_index >= pFinalColorMap->colors)
					hist.histIterms[j].tmp.likely_colormap_index = 0; // actual value doesn't matter, as the guess is out of date anyway
			}
		}
		palette_error = least_error;
		return 0;
	}

	double createIndexImg(ColorMap* pFinalColorMap, const vector<ARGB>& pixels, const UINT& width, vector<vector_fixed<float, 4> > & finalPalette)
	{
		const float gamma = 1 / 2.2f;
		double remapping_error = 0;

		const UINT nMaxColors = pFinalColorMap->colors;
		auto gamma_lut = make_unique<float[]>(nMaxColors);
		toFloatSetGamma(gamma_lut.get(), gamma, nMaxColors);
		finalPalette.clear();
		finalPalette.resize(nMaxColors);
		auto palette = pFinalColorMap->palette.get();
		for (UINT x = 0; x < pFinalColorMap->colors; ++x) {
			floatToRgba(gamma, palette[x].fcolor, finalPalette[x]);
			rgbaToFloat(gamma_lut.get(), finalPalette[x], palette[x].fcolor); /* saves rounding error introduced by to_rgb, which makes remapping & dithering more accurate */
		}

		NearestMap n;
		nearestInit(*pFinalColorMap, n);
		//const UINT max_threads = omp_get_max_threads();
		auto average_color = make_unique<ViterState[]>(nMaxColors);

		const UINT height = pixels.size() / width;
		UINT pixelIndex = 0;
		for (UINT row = 0; row < height; ++row) {
			UINT last_match = 0;
			for (UINT col = 0; col < width; ++col) {
				FloatPixel px;
				rgbaToFloat(gamma_lut.get(), pixels[pixelIndex++], px);
				float diff = 0;
				last_match = nearestSearch(n, px, last_match, diff);

				remapping_error += diff;
				viterUpdateColor(px, 1.0, *pFinalColorMap, last_match, average_color.get());
			}
		}

		viterFinalize(*pFinalColorMap, average_color.get());
		return remapping_error / pixels.size();
	}

	unsigned short nearestColorIndex(const ARGB* pPalette, const unsigned short nMaxColors, ARGB argb, const UINT pos)
	{
		auto got = nearestMap.find(argb);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= 0)
			c = m_transparentColor;

		double mindist = INT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = 0; i < nMaxColors; ++i) {
			Color c2(pPalette[i]);
			double curdist = sqr(c2.GetA() - c.GetA());
			if (curdist > mindist)
				continue;

			getLab(c2, lab2);
			if (nMaxColors <= 4) {
				curdist += sqr(c2.GetR() - c.GetR());
				if (curdist > mindist)
					continue;

				curdist += sqr(c2.GetG() - c.GetG());
				if (curdist > mindist)
					continue;

				curdist += sqr(c2.GetB() - c.GetB());
				if (hasSemiTransparency) {
					if (curdist > mindist)
						continue;
					curdist += sqr(c2.GetA() - c.GetA());
				}
			}
			else if (nMaxColors > 32 || hasSemiTransparency) {
				if(hasSemiTransparency)
					curdist /= exp(0.75);

				curdist += abs(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += _sqrt(sqr(lab2.A - lab1.A) + sqr(lab2.B - lab1.B));
			}
			else {
				auto deltaL_prime_div_k_L_S_L = CIELABConvertor::L_prime_div_k_L_S_L(lab1, lab2);
				curdist += sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double a1Prime, a2Prime, CPrime1, CPrime2;
				auto deltaC_prime_div_k_L_S_L = CIELABConvertor::C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double barCPrime, barhPrime;
				auto deltaH_prime_div_k_L_S_L = CIELABConvertor::H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist += sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				curdist += CIELABConvertor::R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);				
			}

			if (curdist > mindist)
				continue;

			mindist = curdist;
			k = i;
		}
		nearestMap[argb] = k;
		return k;
	}

	unsigned short closestColorIndex(const ARGB* pPalette, const unsigned short nMaxColors, ARGB argb, const UINT pos)
	{
		UINT k = 0;
		Color c(argb);
		if (c.GetA() <= 0xF)
			c = m_transparentColor;

		vector<unsigned short> closest(5);
		auto got = closestMap.find(argb);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = SHRT_MAX;

			for (; k < nMaxColors; ++k) {
				Color c2(pPalette[k]);
				closest[4] = abs(c.GetA() - c2.GetA()) + abs(c.GetR() - c2.GetR()) + abs(c.GetG() - c2.GetG()) + abs(c.GetB() - c2.GetB());
				if (closest[4] < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = closest[4];
				}
				else if (closest[4] < closest[3]) {
					closest[1] = k;
					closest[3] = closest[4];
				}
			}

			if (closest[3] == SHRT_MAX)
				closest[2] = 0;
		}
		else
			closest = got->second;

		if (closest[2] == 0 || (rand() % (closest[3] + closest[2])) <= closest[3])
			k = closest[0];
		else
			k = closest[1];

		closestMap[argb] = closest;
		return k;
	}

	inline int GetColorIndex(const Color& c)
	{
		return GetARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex >= 0);
	}

	bool quantize_image(const ARGB* pixels, const ColorPalette* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither)
			return dither_image(pixels, pPalette->Entries, nMaxColors, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, qPixels, width, height);

		DitherFn ditherFn = (m_transparentPixelIndex >= 0 || nMaxColors < 256) ? nearestColorIndex : closestColorIndex;
		UINT pixelIndex = 0;
		for (UINT j = 0; j < height; ++j) {
			for (UINT i = 0; i < width; ++i)
				qPixels[pixelIndex++] = ditherFn(pPalette->Entries, nMaxColors, pixels[pixelIndex], i + j);
		}

		BlueNoise::dither(width, height, pixels, pPalette->Entries, nMaxColors, ditherFn, GetColorIndex, qPixels);
		return true;
	}

	int MedianCut::quantizeImg(const vector<ARGB>& pixels, const UINT& width, Mat<float>& saliencyMap_float, ColorPalette* pPalette, UINT& newcolors)
	{
		if (pPalette->Count <= 2) {
			if (m_transparentPixelIndex >= 0) {
				pPalette->Entries[0] = m_transparentColor;
				pPalette->Entries[1] = Color::Black;
			}
			else {
				pPalette->Entries[1] = Color::Black;
				pPalette->Entries[0] = Color::White;
			}
			return 0;
		}

		unique_ptr<ColorMap> pFinalColorMap;
		double palette_error;
		Histogram hist(pixels, width, saliencyMap_float);
		iterMedianCut(pFinalColorMap, hist, newcolors, palette_error);

		const double max_mse = INT_MAX;
		const double iteration_limit = sqr(1 / 1024.0);
		UINT iterations = 17;

		if (!iterations && palette_error < 0 && max_mse < INT_MAX)
			iterations = 1; // otherwise total error is never calculated and MSE limit won't work

		if (iterations) {
			double previous_palette_error = INT_MAX;

			for (UINT i = 0; i < iterations; ++i) {
				palette_error = viterDoIteration(hist, *pFinalColorMap, true);

				if (fabs(previous_palette_error - palette_error) < iteration_limit)
					break;

				if (palette_error > max_mse * 1.5) { // probably hopeless
					if (palette_error > max_mse * 3.0)
						break; // definitely hopeless
					++iterations;
				}

				previous_palette_error = palette_error;
			}
		}

		vector<vector_fixed<float, 4> > finalPalette(pPalette->Count);
		createIndexImg(pFinalColorMap.get(), pixels, width, finalPalette);
		newcolors = finalPalette.size();
		for (UINT k = 0; k < newcolors; ++k) {
			auto pixel = finalPalette[k];
			pPalette->Entries[k] = Color::MakeARGB(static_cast<BYTE>(BYTE_MAX * pixel[3]), static_cast<BYTE>(pixel[0]), static_cast<BYTE>(pixel[1]), static_cast<BYTE>(pixel[2]));

			if (m_transparentPixelIndex >= 0 && pPalette->Entries[k] == m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[k]);
		}

		if (palette_error > max_mse)
			return -1;

		return 0;
	}

	bool MedianCut::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();
		const auto area = (size_t) (bitmapWidth * bitmapHeight);

		vector<ARGB> pixels(area);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor, 0xF, nMaxColors);

		if (nMaxColors > 256)
			nMaxColors = 256;

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors == 256 && pDest->GetPixelFormat() != PixelFormat8bppIndexed)
			pDest->ConvertFormat(PixelFormat8bppIndexed, DitherTypeSolid, PaletteTypeCustom, pPalette, 0);

		if (nMaxColors > 2) {
			UINT pixelIndex = 0;
			Mat<float> saliencyMap(bitmapHeight, bitmapWidth);
			float saliencyBase = 0.1;
			for (UINT y = 0; y < bitmapHeight; ++y) {
				for (UINT x = 0; x < bitmapWidth; ++x) {
					CIELABConvertor::Lab lab1;
					getLab(pixels[pixelIndex++], lab1);
					saliencyMap(y, x) = saliencyBase + (1 - saliencyBase) * lab1.L / 255.0f;
				}
			}
			quantizeImg(pixels, bitmapWidth, saliencyMap, pPalette, nMaxColors);
		}
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette->Entries[0] = m_transparentColor;
				pPalette->Entries[1] = Color::Black;
			}
			else {
				pPalette->Entries[0] = Color::Black;
				pPalette->Entries[1] = Color::White;
			}
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		if (hasSemiTransparency || nMaxColors <= 32)
			PR = PG = PB = 1;
		quantize_image(pixels.data(), pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);

		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if (nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}
		pixelMap.clear();
		closestMap.clear();
		nearestMap.clear();

		pDest->SetPalette(pPalette);
		return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex >= 0);
	}
}
