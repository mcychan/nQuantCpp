#pragma once
/* Copyright (c) 2006 Derrick Coetzee
Copyright(c) 2015 Hao-Zhi Huang
Copyright (c) 2018 Miller Cy Chan

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "stdafx.h"
#include "EdgeAwareSQuantizer.h"
#include "PnnLABQuantizer.h"

#include <deque>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <time.h>
#include <limits>
#include <map>

namespace EdgeAwareSQuant
{
	bool hasTransparency = false;
	ARGB m_transparentColor = Color::Transparent;

	const int DECOMP_SVD = 1;

	static bool mycmp(pair<float, int> p1, pair<float, int> p2)
	{
		return p1.first < p2.first;
	}

	template <typename T>
	class array3d
	{
	public:
		array3d(int width, int height, int depth)
		{
			this->width = width;
			this->height = height;
			this->depth = depth;
			data = make_unique<T[]>(width * height * depth);
		}

		array3d(const array3d<T>& rhs)
		{
			width = rhs.get_width();
			height = rhs.height;
			depth = rhs.depth;
			data = make_unique<T[]>(width * height * depth);
			for (int i = 0; i< (width * height * depth); i++)
				data[i] = rhs.data[i];
		}

		inline T& operator()(int col, int row, int layer)
		{
			return data[row * width * depth + col * depth + layer];
		}

		inline int get_width()  const { return width; }
		inline int get_height() const { return height; }
		inline int get_depth() const { return depth; }

	private:
		unique_ptr<T[]> data;
		int width, height, depth;
	};

	int compute_max_coarse_level(int width, int height) {
		// We want the coarsest layer to have at most MAX_PIXELS pixels
		const int MAX_PIXELS = 4000;
		int result = 0;
		while (width * height > MAX_PIXELS) {
			width >>= 1;
			height >>= 1;
			result++;
		}
		return result;
	}

	void fill_random_icm(Mat<byte>& indexImg8, int palette_size) {
		for (int i = 0; i < indexImg8.get_height(); i++) {
			for (int j = 0; j < indexImg8.get_width(); j++) {
				int ran_val = float(rand()) / RAND_MAX * (palette_size - 1);
				if (ran_val < 0)
					ran_val = 0;
				if (ran_val >= palette_size)
					ran_val = palette_size - 1;
				indexImg8(i, j) = ran_val;
			}
		}
	}

	void random_permutation(int count, vector<int>& result) {
		result.resize(count);
		iota(result.begin(), result.end(), 0);
		random_shuffle(result.begin(), result.end());
	}

	void random_permutation_2d(int width, int height, deque<pair<int, int> >& result) {
		vector<int> perm1d;
		random_permutation(width * height, perm1d);
		for (auto it = perm1d.rbegin(); it != perm1d.rend(); ++it)
			result.emplace_back(*it % width, *it / width);
	}

	void compute_b_array_ea_saliency(Mat<Mat<float> >& weightMaps, Mat<Mat<float> >& b, int filterRadius, Mat<float>& saliencyMap)
	{
		int imgHeight = weightMaps.get_height();
		int imgWidth = weightMaps.get_width();
		for (int i_y = 0; i_y < imgHeight; i_y++) {
			for (int i_x = 0; i_x < imgWidth; i_x++) {
				auto& wmap_i = weightMaps(i_y, i_x);
				int extendedFilterRadius = filterRadius * 2;
				auto& b_yx = b(i_y, i_x);
				b_yx.reset(extendedFilterRadius * 2 + 1, extendedFilterRadius * 2 + 1);
				int j_y_min = i_y - extendedFilterRadius, j_y_max = i_y + extendedFilterRadius;
				int j_x_min = i_x - extendedFilterRadius, j_x_max = i_x + extendedFilterRadius;
				int wmI_y_min = i_y - filterRadius, wmI_y_max = i_y + filterRadius, wmI_x_min = i_x - filterRadius, wmI_x_max = i_x + filterRadius;
				for (int j_y = j_y_min; j_y <= j_y_max; j_y++) {
					if (j_y < 0 || j_y >= imgHeight)
						continue;
					for (int j_x = j_x_min; j_x <= j_x_max; j_x++) {
						if (j_x < 0 || j_x >= imgWidth)
							continue;

						auto& wmap_j = weightMaps(j_y, j_x);
						int wmJ_y_min = j_y - filterRadius, wmJ_y_max = j_y + filterRadius, wmJ_x_min = j_x - filterRadius, wmJ_x_max = j_x + filterRadius;
						for (int wmJ_y = wmJ_y_min; wmJ_y <= wmJ_y_max; wmJ_y++) {
							if (wmJ_y < 0 || wmJ_y >= imgHeight)
								continue;
							for (int wmJ_x = wmJ_x_min; wmJ_x <= wmJ_x_max; wmJ_x++) {
								if (wmJ_x < 0 || wmJ_x >= imgWidth)
									continue;
								// if in overlap area
								if (abs(wmJ_y - i_y) <= filterRadius && abs(wmJ_x - i_x) <= filterRadius)
									b_yx(j_y - j_y_min, j_x - j_x_min) += saliencyMap(wmJ_y, wmJ_x) * wmap_i(wmJ_y - wmI_y_min, wmJ_x - wmI_x_min) * wmap_j(wmJ_y - wmJ_y_min, wmJ_x - wmJ_x_min);
							}
						}

						if (b_yx(j_y - j_y_min, j_x - j_x_min) == 0)
							b_yx(j_y - j_y_min, j_x - j_x_min) = 1e-10;
					}
				}
			}
		}
	}

	float b_value_ea(Mat<Mat<float> >& b, int i_x, int i_y, int j_x, int j_y)
	{
		auto& b_yx = b(i_y, i_x);
		int extendedFilterRadius = (b_yx.get_height() - 1) / 2;
		int k_x = j_x - i_x + extendedFilterRadius;
		int k_y = j_y - i_y + extendedFilterRadius;
		if (k_x < 0 || k_y < 0 || k_x >= b_yx.get_width() || k_y >= b_yx.get_height())
			return 1e-10;
		return b_yx(k_y, k_x);
	}

	void compute_a_image_ea(vector<ARGB>& image, Mat<Mat<float> >& b, array2d<vector_fixed<float, 3> >& a)
	{
		int extendedFilterRadius = (b(0, 0).get_width() - 1) / 2;
		for (int i_y = 0; i_y < a.get_height(); i_y++) {
			for (int i_x = 0; i_x < a.get_width(); i_x++) {
				for (int j_y = i_y - extendedFilterRadius; j_y <= i_y + extendedFilterRadius; j_y++) {
					if (j_y < 0)
						j_y = 0;
					if (j_y >= a.get_height())
						break;

					for (int j_x = i_x - extendedFilterRadius; j_x <= i_x + extendedFilterRadius; j_x++) {
						if (j_x < 0)
							j_x = 0;
						if (j_x >= a.get_width())
							break;

						float tmpBvalue = b_value_ea(b, i_x, i_y, j_x, j_y);
						if (tmpBvalue != 0) {
							Color jPixel(image[j_y * a.get_width() + j_x]);
							a(i_x, i_y)(0) += tmpBvalue * jPixel.GetR() / 255.0;
							a(i_x, i_y)(1) += tmpBvalue * jPixel.GetG() / 255.0;
							a(i_x, i_y)(2) += tmpBvalue * jPixel.GetB() / 255.0;
						}
					}
				}
				a(i_x, i_y) *= -2.0;
			}
		}
	}

	template <typename T, int length>
	array2d<T> extract_vector_layer_2d(array2d<vector_fixed<T, length> >& s, int k)
	{
		array2d<T> result(s.get_width(), s.get_height());
		for (int i = 0; i < s.get_width(); i++) {
			for (int j = 0; j < s.get_height(); j++)
				result(i, j) = s(i, j)(k);
		}
		return result;
	}

	template <typename T, int length>
	vector<T> extract_vector_layer_1d(vector<vector_fixed<T, length> >& s, int k)
	{
		vector<T> result(s.size());
		for (UINT i = 0; i < s.size(); i++)
			result[i] = s[i](k);

		return result;
	}

	void sum_coarsen(array2d<vector_fixed<float, 3> >& fine, array2d<vector_fixed<float, 3> >& coarse)
	{
		for (int y = 0; y<coarse.get_height(); y++) {
			for (int x = 0; x<coarse.get_width(); x++) {
				coarse(x, y) = fine(x * 2, y * 2);
				if (x * 2 + 1 < fine.get_width())
					coarse(x, y) += fine(x * 2 + 1, y * 2);
				if (y * 2 + 1 < fine.get_height())
					coarse(x, y) += fine(x * 2, y * 2 + 1);
				if (x * 2 + 1 < fine.get_width() && y * 2 + 1 < fine.get_height())
					coarse(x, y) += fine(x * 2 + 1, y * 2 + 1);
			}
		}
	}

	void zoom_float_icm(Mat<byte>& smallVal, Mat<byte>& big)
	{
		for (int y = 0; y < big.get_height(); y++) {
			int small_y = y / 2.0;
			if (small_y >= smallVal.get_height())
				continue;
			for (int x = 0; x < big.get_width(); x++) {
				int small_x = x / 2.0;
				if (small_x >= smallVal.get_width())
					continue;
				big(y, x) = smallVal(small_y, small_x);
			}
		}
	}

	void compute_initial_s_ea_icm(array2d<vector_fixed<float, 3> >& s, Mat<byte>& indexImg8, Mat<Mat<float> >& b)
	{
		int palette_size = s.get_width();
		int coarse_width = indexImg8.get_width();
		int coarse_height = indexImg8.get_height();
		int center_x = (b(0, 0).get_height() - 1) / 2, center_y = (b(0, 0).get_height() - 1) / 2;
		int extendedFilterRadius = (b(0, 0).get_height() - 1) / 2;
		vector_fixed<float, 3> zero_vector;
		for (int v = 0; v < palette_size; v++) {
			for (int alpha = v; alpha < palette_size; alpha++)
				s(v, alpha) = zero_vector; // alpha > v
		}
		for (int i_y = 0; i_y < coarse_height; i_y++) {
			for (int i_x = 0; i_x < coarse_width; i_x++) {
				int j_y_min = i_y - extendedFilterRadius, j_y_max = i_y + extendedFilterRadius;
				int j_x_min = i_x - extendedFilterRadius, j_x_max = i_x + extendedFilterRadius;
				for (int j_y = j_y_min; j_y <= j_y_max; j_y++) {
					if (j_y < 0 || j_y >= coarse_height)
						continue;
					for (int j_x = j_x_min; j_x <= j_x_max; j_x++) {
						if (j_x < 0 || j_x >= coarse_width)
							continue;
						if (i_x == j_x && i_y == j_y)
							continue;
						float b_ij = b_value_ea(b, i_x, i_y, j_x, j_y);
						int v = indexImg8(i_y, i_x);
						int alpha = indexImg8(j_y, j_x);
						s(v, alpha)(0) += b_ij;
						s(v, alpha)(1) += b_ij;
						s(v, alpha)(2) += b_ij;
					}
				}
				int v = indexImg8(i_y, i_x);
				float b_ii = b_value_ea(b, i_x, i_y, i_x, i_y);
				s(v, v)(0) += b_ii;
				s(v, v)(1) += b_ii;
				s(v, v)(2) += b_ii;
			}
		}
	}

	void refine_palette_icm_mat(array2d<vector_fixed<float, 3> >& s, Mat<byte>& indexImg8,
		array2d<vector_fixed<float, 3> >& a, vector<vector_fixed<float, 3> >& palette, int& palatte_changed)
	{
		// We only computed the half of S above the diagonal - reflect it
		for (int v = 0; v < s.get_width(); v++) {
			for (int alpha = 0; alpha < v; alpha++) {
				s(v, alpha) = s(alpha, v);
				ASSERT(!isnan(s(v, alpha)(0)));
				ASSERT(!isnan(s(v, alpha)(1)));
				ASSERT(!isnan(s(v, alpha)(2)));
			}
			ASSERT(!isnan(s(v, v)(0)));
			ASSERT(!isnan(s(v, v)(1)));
			ASSERT(!isnan(s(v, v)(2)));
		}

		vector<vector_fixed<float, 3> > r(palette.size());

		for (int i_y = 0; i_y < indexImg8.get_height(); i_y++) {
			for (int i_x = 0; i_x < indexImg8.get_width(); i_x++) {
				int label = indexImg8(i_y, i_x);
				r[label] += a(i_x, i_y);
			}
		}

		float max_palette_delta = 0, min_palette_delta = 1.0;
		for (UINT k = 0; k < 3; k++) {
			bool isBlack = true;
			auto& S_k = extract_vector_layer_2d(s, k);
			auto& R_k = extract_vector_layer_1d(r, k);
			auto& palette_channel = -1.0f * ((2.0f * S_k).matrix_inverse()) * R_k;
			UINT v = 0;
			for (; v < palette.size(); v++) {
				float val = palette_channel[v];
				if (val < 0.0 || isnan(val))
					val = 0.0;
				if (val > 1.0)
					val = 1.0;
				if (val != 0)
					isBlack = false;
				float palette_delta = abs(palette[v](k) - val);
				if (palette_delta > max_palette_delta)
					max_palette_delta = palette_delta;
				if (palette_delta > min_palette_delta)
					min_palette_delta = palette_delta;
				if (palette_delta > 1.0 / 255.0)
					palatte_changed++;
				palette[v](k) = val;
			}
			if (isBlack)
				swap(palette[v], palette[0]);
		}
	}

	void spatial_color_quant_ea_icm_saliency(vector<ARGB>& image, Mat<Mat<float> >& weightMaps, Mat<float> saliencyMap,
		array2d<UINT>& quantized_image, vector<vector_fixed<float, 3> >& palette,
		const float initial_temperature = 1.0, const float final_temperature = 0.00001, const int temps_per_level = 1, const int repeats_per_temp = 1, const int filter_radius = 1)
	{
		int allNeiLevel = 1;
		int max_coarse_level = 4;
		int neiSize = 10;

		auto ppIndexImg8 = make_unique<Mat<byte> >(weightMaps.get_height() >> max_coarse_level, weightMaps.get_width() >> max_coarse_level);
		auto pIndexImg8 = ppIndexImg8.get();
		fill_random_icm(*pIndexImg8, palette.size());

		// Compute a_I^l, b_{IJ}^l according to  Puzicha's (18)
		auto a_array = make_unique<array2d<vector_fixed<float, 3> >[]>(max_coarse_level + 1);
		auto b_array = make_unique<Mat<Mat<float> >[]>(max_coarse_level + 1);

		auto& b0 = b_array[0];
		b0.reset(weightMaps.get_height(), weightMaps.get_width());
		compute_b_array_ea_saliency(weightMaps, b0, filter_radius, saliencyMap);

		auto& a0 = a_array[0];
		a0.reset(weightMaps.get_width(), weightMaps.get_height());
		compute_a_image_ea(image, b0, a0);

		int coarse_level = 1;
		for (; coarse_level <= max_coarse_level; coarse_level++) {
			auto& ai = a_array[coarse_level];
			ai.reset(weightMaps.get_width() >> coarse_level, weightMaps.get_height() >> coarse_level);

			int newExtendedFilterSize = b0(0, 0).get_height() - 2;
			newExtendedFilterSize = max(3, newExtendedFilterSize);

			auto& bi = b_array[coarse_level];
			bi.reset(ai.get_height(), ai.get_width());
			int newExtendedFilterRadius = (newExtendedFilterSize - 1) / 2;

			for (int I_y = 0; I_y < ai.get_height(); I_y++) {
				for (int I_x = 0; I_x < ai.get_width(); I_x++) {
					auto& bi_yx = bi(I_y, I_x);
					bi_yx.reset(newExtendedFilterSize, newExtendedFilterSize);
					int J_y_min = I_y - newExtendedFilterRadius, J_y_max = I_y + newExtendedFilterRadius;
					int J_x_min = I_x - newExtendedFilterRadius, J_x_max = I_x + newExtendedFilterRadius;

					for (int J_y = J_y_min; J_y <= J_y_max; J_y++) {
						if (J_y < 0 || J_y >= ai.get_height())
							continue;
						for (int J_x = J_x_min; J_x <= J_x_max; J_x++) {
							if (J_x < 0 || J_x >= ai.get_width())
								continue;

							for (int i_y = I_y * 2; i_y < I_y * 2 + 2; i_y++) {
								if (i_y >= a0.get_height())
									continue;
								for (int i_x = I_x * 2; i_x < I_x * 2 + 2; i_x++) {
									if (i_x >= a0.get_width())
										continue;
									for (int j_y = J_y * 2; j_y < J_y * 2 + 2; j_y++) {
										if (j_y >= a0.get_height())
											continue;
										for (int j_x = J_x * 2; j_x < J_x * 2 + 2; j_x++) {
											if (j_x >= a0.get_width())
												continue;

											bi_yx(J_y - J_y_min, J_x - J_x_min) += b_value_ea(b0, i_x, i_y, j_x, j_y);
										}
									}
								}
							}

							if (bi_yx(J_y - J_y_min, J_x - J_x_min) == 0)
								bi_yx(J_y - J_y_min, J_x - J_x_min) = 1e-10;

						}
					}
				}
			}

			sum_coarsen(a0, ai);
		}


		// Multiscale ICM
		coarse_level = max_coarse_level;
		array2d<vector_fixed<float, 3> > s(palette.size(), palette.size());
		compute_initial_s_ea_icm(s, *pIndexImg8, b_array[coarse_level]);

		float paletteSize = palette.size() * 1.0;
		while (coarse_level >= 0) {
			// calculate the distance between centroids
			vector<vector<pair<float, int> > > centroidDist(paletteSize, vector<pair<float, int> >(paletteSize, pair<float, int>(0.0f, -1)));
			for (int l1 = 0; l1 < palette.size(); l1++) {
				for (int l2 = l1; l2 < palette.size(); l2++) {
					float curDist = 0.0;
					for (int c = 0; c < 3; c++)
						curDist += pow(palette[l1](c) - palette[l2](c), 2);

					centroidDist[l1][l2] = pair<float, int>(curDist, l2);
					centroidDist[l2][l1] = pair<float, int>(curDist, l1);
				}
			}
			// sort centroidDist row by row
			for (int l1 = 0; l1 < palette.size(); l1++)
				sort(centroidDist[l1].begin(), centroidDist[l1].end(), mycmp);

			auto& a = a_array[coarse_level];
			auto& b = b_array[coarse_level];
			const auto& b1 = b(0, 0);

			int center_x = (b1.get_width() - 1) / 2, center_y = (b1.get_height() - 1) / 2;

			int step_counter = 0;
			int repeat_outter = 0;
			int palette_changed = 0;
			while (repeat_outter == 0 || palette_changed > palette.size() * 0.1) {
				palette_changed = 0;
				repeat_outter++;
				//----update labeling
				int pixels_changed = 0, pixels_visited = 0;
				int repeat_inner = 0;

				while (repeat_inner == 0 || pixels_changed > 0.001 * pIndexImg8->get_width() * pIndexImg8->get_height()) {
					repeat_inner++;
					pixels_changed = 0;
					pixels_visited = 0;

					deque<pair<int, int> > visit_queue;
					random_permutation_2d(pIndexImg8->get_width(), pIndexImg8->get_height(), visit_queue);

					// Compute 2*sum(j in extended neighborhood of i, j != i) b_ij
					while (visit_queue.size() > 0.0 * pIndexImg8->get_width() * pIndexImg8->get_height()) {
						// pick a pixel every time
						int i_x = visit_queue.front().first, i_y = visit_queue.front().second;
						visit_queue.pop_front();

						// Compute based on Puzicha's (28)
						vector_fixed<float, 3> p_i;
						for (int j_y = i_y - (b1.get_height() - 1) / 2; j_y <= i_y + (b1.get_height() - 1) / 2; j_y++) {
							if (j_y < 0 || j_y >= pIndexImg8->get_height())
								continue;
							for (int j_x = i_x - (b1.get_width() - 1) / 2; j_x <= i_x + (b1.get_width() - 1) / 2; j_x++) {
								//int j_x = x - center_x + i_x, j_y = y - center_y + i_y;
								if (i_x == j_x && i_y == j_y)
									continue;
								if (j_x < 0 || j_x >= pIndexImg8->get_width())
									continue;
								float b_ij = b_value_ea(b, i_x, i_y, j_x, j_y);
								auto& pixelIndex = pIndexImg8->at(j_y, j_x);
								p_i(0) += b_ij * palette[pixelIndex](0);
								p_i(1) += b_ij * palette[pixelIndex](1);
								p_i(2) += b_ij * palette[pixelIndex](2);
							}
						}

						p_i *= 2.0;
						p_i += a(i_x, i_y);

						int old_max_v = pIndexImg8->at(i_y, i_x);

						vector<float> meanfields;
						float min_meanfield = (numeric_limits<float>::max)();
						float middle_b = b_value_ea(b, i_x, i_y, i_x, i_y);
						int bestLabel = old_max_v;

						if (coarse_level >= allNeiLevel) {
							// search for all palette color
							for (UINT v = 0; v < palette.size(); v++) {
								float mf_val = palette[v].dot_product(p_i + palette[v] * middle_b);
								if (mf_val < min_meanfield) {
									min_meanfield = mf_val;
									bestLabel = v;
								}
							}
						}
						else {
							// just looking for the palette color which is near current palette color
							if (neiSize > palette.size())
								neiSize = palette.size();
							for (int v = 0; v < neiSize; ++v) {
								int tryLabel = centroidDist[old_max_v][v].second;
								float mf_val = palette[tryLabel].dot_product(p_i + palette[tryLabel] * middle_b);
								if (mf_val < min_meanfield) {
									min_meanfield = mf_val;
									bestLabel = tryLabel;
								}
							}
						}


						pIndexImg8->at(i_y, i_x) = bestLabel;
						if ((palette[bestLabel] - palette[old_max_v]).norm_squared() >= 1.0 / (255.0 * 255.0))
							pixels_changed++;

						pixels_visited++;
					}
				}

				//----update palette----
				compute_initial_s_ea_icm(s, *pIndexImg8, b_array[coarse_level]);
				refine_palette_icm_mat(s, *pIndexImg8, a, palette, palette_changed);
			}

			if (--coarse_level < 0)
				break;
			auto pNewIndexImg8 = new Mat<byte>(weightMaps.get_height() >> coarse_level, weightMaps.get_width() >> coarse_level);
			zoom_float_icm(*pIndexImg8, *pNewIndexImg8);
			ppIndexImg8.reset(pNewIndexImg8);
			pIndexImg8 = ppIndexImg8.get();
		}

		a_array.reset();
		b_array.reset();

		for (int i_x = 0; i_x < weightMaps.get_width(); i_x++) {
			for (int i_y = 0; i_y < weightMaps.get_height(); i_y++)
				quantized_image(i_x, i_y) = pIndexImg8->at(i_y, i_x);
		}
	}

	void filter_bila(vector<ARGB>& img, Mat<Mat<float> >& weightMaps, const float sigma_s = 1.0f, const float sigma_r = 2.0f, const int iteration = 4)
	{
		// pixel-wise filter		
		int radius = 1;
		float wMin = 100;
		for (int y = 0; y < weightMaps.get_height(); y++) {
			for (int x = 0; x < weightMaps.get_width(); x++) {
				float weightSum = 0.0;

				float sum[3] = { 0 };
				int yyMin = y - radius, yyMax = y + radius, xxMin = x - radius, xxMax = x + radius;
				auto& weightMaps_yx = weightMaps(y, x);
				weightMaps_yx.reset(2 * radius + 1, 2 * radius + 1);
				for (int yy = yyMin; yy <= yyMax; yy++) {
					if (yy < 0 || yy >= weightMaps.get_height())
						continue;
					for (int xx = xxMin; xx <= xxMax; xx++) {
						if (xx < 0 || xx >= weightMaps.get_width())
							continue;
						float spaceD = pow(y - yy, 2) + pow(x - xx, 2);
						Color pixelXY(img[y * weightMaps.get_width() + x]);
						Color pixelXXYY(img[yy * weightMaps.get_width() + xx]);
						float colorD = pow(pixelXY.GetR() - pixelXXYY.GetR(), 2) + pow(pixelXY.GetG() - pixelXXYY.GetG(), 2) + pow(pixelXY.GetB() - pixelXXYY.GetB(), 2);
						float tmpW = BYTE_MAX * exp(-spaceD / (2 * sigma_s * sigma_s) - colorD / (2 * sigma_r * sigma_r));

						weightSum += tmpW;

						sum[0] += tmpW * pixelXXYY.GetR() / 255.0;
						sum[1] += tmpW * pixelXXYY.GetG() / 255.0;
						sum[2] += tmpW * pixelXXYY.GetB() / 255.0;
						weightMaps_yx(xx - xxMin, yy - yyMin) = tmpW;

						if (tmpW < wMin)
							wMin = tmpW;
					}
				}

				weightMaps_yx /= weightSum;
			}
		}
	}

	bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, array2d<UINT>& quantized_image)
	{
		pDest->SetPalette(pPalette);

		BitmapData targetData;
		UINT w = pDest->GetWidth();
		UINT h = pDest->GetHeight();

		Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
		if (status != Ok) {
			AfxMessageBox(_T("Cannot write image"));
			return false;
		}

		auto pRowDest = (byte*)targetData.Scan0;
		UINT strideDest;

		// Compensate for possible negative stride
		if (targetData.Stride > 0)
			strideDest = targetData.Stride;
		else {
			pRowDest += h * targetData.Stride;
			strideDest = -targetData.Stride;
		}

		UINT bpp = GetPixelFormatSize(pDest->GetPixelFormat());
		// Second loop: fill indexed bitmap
		for (UINT y = 0; y < h; y++) {	// For each row...
			for (UINT x = 0; x < w; x++) {	// ...for each pixel...
				byte nibbles = 0;
				byte index = static_cast<byte>(quantized_image(x, y));

				switch (bpp)
				{
				case 8:
					pRowDest[x] = index;
					break;
				case 4:
					// First pixel is the high nibble. From and To indices are 0..16
					nibbles = pRowDest[x / 2];
					if ((x & 1) == 0) {
						nibbles &= 0x0F;
						nibbles |= (byte)(index << 4);
					}
					else {
						nibbles &= 0xF0;
						nibbles |= index;
					}

					pRowDest[x / 2] = nibbles;
					break;
				case 1:
					// First pixel is MSB. From and To are 0 or 1.
					int pos = x / 8;
					byte mask = (byte)(128 >> (x & 7));
					if (index == 0)
						pRowDest[pos] &= (byte)~mask;
					else
						pRowDest[pos] |= mask;
					break;
				}
			}

			pRowDest += strideDest;
		}

		status = pDest->UnlockBits(&targetData);
		return pDest->GetLastStatus() == Ok;
	}

	bool EdgeAwareSQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		if (nMaxColors > 256)
			nMaxColors = 256;
		UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		UINT bitmapWidth = pSource->GetWidth();
		UINT bitmapHeight = pSource->GetHeight();

		hasTransparency = false;
		bool r = true;
		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++) {
					Color color;
					pSource->GetPixel(x, y, &color);

					if (color.GetA() < BYTE_MAX) {
						hasTransparency = true;
						if (color.GetA() == 0)
							m_transparentColor = color.GetValue();
					}
					pixels[pixelIndex++] = color.GetValue();
				}
			}
		}

		// Lock bits on 3x8 source bitmap
		else {
			BitmapData data;
			Status status = pSource->LockBits(&Rect(0, 0, bitmapWidth, bitmapHeight), ImageLockModeRead, pSource->GetPixelFormat(), &data);
			if (status != Ok)
				return false;

			auto pRowSource = (byte*)data.Scan0;
			UINT strideSource;

			if (data.Stride > 0) strideSource = data.Stride;
			else
			{
				// Compensate for possible negative stride
				// (not needed for first loop, but we have to do it
				// for second loop anyway)
				pRowSource += bitmapHeight * data.Stride;
				strideSource = -data.Stride;
			}

			int pixelIndex = 0;

			// First loop: gather color information
			for (UINT y = 0; y < bitmapHeight; y++) {	// For each row...
				auto pPixelSource = pRowSource;

				for (UINT x = 0; x < bitmapWidth; x++) {	// ...for each pixel...
					byte pixelBlue = *pPixelSource++;
					byte pixelGreen = *pPixelSource++;
					byte pixelRed = *pPixelSource++;
					byte pixelAlpha = bitDepth < 32 ? BYTE_MAX : *pPixelSource++;

					auto argb = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
					if (pixelAlpha < BYTE_MAX) {
						hasTransparency = true;
						if (pixelAlpha == 0)
							m_transparentColor = argb;
					}
					pixels[pixelIndex++] = argb;
				}

				pRowSource += strideSource;
			}

			pSource->UnlockBits(&data);
		}

		// see equation (7) in the paper
		Mat<float> saliencyMap(bitmapHeight, bitmapWidth);
		float saliencyBase = 0.1;
		for (int y = 0; y < saliencyMap.get_height(); y++) {
			for (int x = 0; x < saliencyMap.get_width(); x++)
				saliencyMap(y, x) = saliencyBase + (1 - saliencyBase) / 255.0;
		}

		array2d<UINT> quantized_image(bitmapWidth, bitmapHeight);
		auto bins = make_unique<PnnLABQuant::pnnbin[]>(65536);
		auto pPaletteBytes = make_unique<byte[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
		pnnLABQuantizer.pnnquan(pixels, bins.get(), pPalette);

		// init
		vector<vector_fixed<float, 3> > palette(nMaxColors);
		for (int k = 0; k < nMaxColors; k++) {
			Color c(pPalette->Entries[k]);
			palette[k](0) = c.GetR() / 255.0;
			palette[k](1) = c.GetG() / 255.0;
			palette[k](2) = c.GetB() / 255.0;
		}

		Mat<Mat<float> > weightMaps(bitmapHeight, bitmapWidth);
		filter_bila(pixels, weightMaps);
		spatial_color_quant_ea_icm_saliency(pixels, weightMaps, saliencyMap, quantized_image, palette);

		if (nMaxColors > 2) {
			/* Fill palette */
			for (UINT k = 0; k<nMaxColors; ++k) {
				pPalette->Entries[k] = Color::MakeARGB(BYTE_MAX, static_cast<byte>(BYTE_MAX * palette[k](0)), static_cast<byte>(BYTE_MAX * palette[k](1)), static_cast<byte>(BYTE_MAX * palette[k](2)));
				if (hasTransparency && pPalette->Entries[k] == m_transparentColor)
					swap(pPalette->Entries[0], pPalette->Entries[k]);
			}
		}
		else {
			if (hasTransparency) {
				pPalette->Entries[1] = Color::Transparent;
				pPalette->Entries[0] = Color::Black;
			}
			else {
				pPalette->Entries[1] = Color::White;
				pPalette->Entries[0] = Color::Black;
			}
		}

		return ProcessImagePixels(pDest, pPalette, quantized_image);
	}

}