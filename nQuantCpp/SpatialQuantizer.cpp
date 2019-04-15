#pragma once
/* Copyright (c) 2006 Derrick Coetzee
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
#include "SpatialQuantizer.h"
#include "bitmapUtilities.h"

#include <deque>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <time.h>
#include <limits>

namespace SpatialQuant
{
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;

	template <typename T, int length>
	class vector_fixed
	{
	public:
		vector_fixed()
		{
			for (int i = 0; i<length; i++)
				data[i] = 0;
		}

		vector_fixed(const vector_fixed<T, length>& rhs)
		{
			for (int i = 0; i<length; i++)
				data[i] = rhs.data[i];
		}

		vector_fixed(const vector<T>& rhs)
		{
			for (int i = 0; i<length; i++)
				data[i] = rhs.data[i];
		}

		inline T& operator[](int index)
		{
			return data[index];
		}

		inline const T& operator[](int index) const
		{
			return data[index];
		}

		inline int get_length() const { return length; }

		T norm_squared() {
			T result = 0;
			for (int i = 0; i<length; i++)
				result += data[i] * data[i];

			return result;
		}

		vector_fixed<T, length>& operator=(const vector_fixed<T, length>& rhs)
		{
			for (int i = 0; i<length; i++)
				data[i] = rhs.data[i];

			return *this;
		}

		vector_fixed<T, length> direct_product(const vector_fixed<T, length>& rhs) {
			vector_fixed<T, length> result;
			for (int i = 0; i<length; i++)
				result[i] = data[i] * rhs.data[i];

			return result;
		}

		T dot_product(const vector_fixed<T, length>& rhs) {
			T result = 0;
			for (int i = 0; i<length; i++) {
				result += data[i] * rhs.data[i];
			}
			return result;
		}

		vector_fixed<T, length>& operator+=(const vector_fixed<T, length>& rhs) {
			for (int i = 0; i<length; i++) {
				data[i] += rhs.data[i];
			}
			return *this;
		}

		vector_fixed<T, length> operator+(const vector_fixed<T, length>& rhs) {
			vector_fixed<T, length> result(*this);
			result += rhs;
			return result;
		}

		vector_fixed<T, length>& operator-=(const vector_fixed<T, length>& rhs) {
			for (int i = 0; i<length; i++) {
				data[i] -= rhs.data[i];
			}
			return *this;
		}

		vector_fixed<T, length> operator-(const vector_fixed<T, length>& rhs) {
			vector_fixed<T, length> result(*this);
			result -= rhs;
			return result;
		}

		vector_fixed<T, length>& operator*=(const T scalar) {
			for (int i = 0; i<length; i++)
				data[i] *= scalar;

			return *this;
		}

		vector_fixed<T, length> operator*(const T scalar) {
			vector_fixed<T, length> result(*this);
			result *= scalar;
			return result;
		}

	private:
		T data[length];
	};

	template <typename T, int length>
	vector_fixed<T, length> operator*(const T scalar, vector_fixed<T, length> vec) {
		return vec * scalar;
	}


	template <typename T>
	class array2d
	{
	public:
		array2d(int width, int height)
		{
			this->width = width;
			this->height = height;
			data = make_unique<T[]>(width * height);
		}

		array2d(const array2d<T>& rhs)
		{
			width = rhs.width;
			height = rhs.height;
			data = make_unique<T[]>(width * height);
			for (int i = 0; i<(width * height); i++)
				data[i] = rhs.data[i];
		}

		inline T& operator()(int col, int row)
		{
			return data[row * width + col];
		}

		inline const T& operator()(int col, int row) const
		{
			return data[row * width + col];
		}

		inline const T& operator[](int index) const
		{
			return data[index];
		}

		inline int get_width() const { return width; }
		inline int get_height() const { return height; }

		array2d<T>& operator*=(const T scalar) {
			for (int i = 0; i<(width * height); i++)
				data[i] *= scalar;
			return *this;
		}

		array2d<T> operator*(const T scalar) {
			array2d<T> result(*this);
			result *= scalar;
			return result;
		}

		vector<T> operator*(const vector<T>& vec) {
			vector<T> result(get_height());
			for (int row = 0; row<get_height(); row++) {
				T sum = 0;
				for (int col = 0; col<get_width(); col++)
					sum += (*this)(col, row) * vec[col];

				result[row] = sum;
			}
			return result;
		}

		array2d<T>& multiply_row_scalar(int row, T mult) {
			for (int i = 0; i<get_width(); i++)
				(*this)(i, row) *= mult;

			return *this;
		}

		array2d<T>& add_row_multiple(int from_row, int to_row, T mult) {
			for (int i = 0; i<get_width(); i++)
				(*this)(i, to_row) += mult * (*this)(i, from_row);

			return *this;
		}

		// We use simple Gaussian elimination - perf doesn't matter since
		// the matrices will be K x K, where K = number of palette entries.
		array2d<T> matrix_inverse() {
			array2d<T> result(get_width(), get_height());
			auto& a = *this;

			// Set result to identity matrix
			for (int i = 0; i<get_width(); i++)
				result(i, i) = 1;

			// Reduce to echelon form, mirroring in result
			for (int i = 0; i<get_width(); i++) {
				result.multiply_row_scalar(i, 1 / a(i, i));
				multiply_row_scalar(i, 1 / a(i, i));
				for (int j = i + 1; j<get_height(); j++) {
					result.add_row_multiple(i, j, -a(i, j));
					add_row_multiple(i, j, -a(i, j));
				}
			}
			// Back substitute, mirroring in result
			for (int i = get_width() - 1; i >= 0; i--) {
				for (int j = i - 1; j >= 0; j--) {
					result.add_row_multiple(i, j, -a(i, j));
					add_row_multiple(i, j, -a(i, j));
				}
			}
			// result is now the inverse
			return result;
		}

	private:
		unique_ptr<T[]> data;
		int width, height;
	};

	template <typename T>
	array2d<T> operator*(T scalar, array2d<T> a) {
		return a * scalar;
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
			width = rhs.width;
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

		inline const T& operator()(int col, int row, int layer) const
		{
			return data[row * width * depth + col * depth + layer];
		}

		void fill_random() {
			for (int i = 0; i< (width * height * depth); i++)
				data[i] = ((double)rand()) / RAND_MAX;
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

	void compute_b_array(array2d<vector_fixed<double, 4> >& filter_weights,
		array2d<vector_fixed<double, 4> >& b)
	{
		// Assume that the pixel i is always located at the center of b,
		// and vary pixel j's location through each location in b.
		int radius_width = (filter_weights.get_width() - 1) / 2,
			radius_height = (filter_weights.get_height() - 1) / 2;
		int offset_x = (b.get_width() - 1) / 2 - radius_width;
		int offset_y = (b.get_height() - 1) / 2 - radius_height;
		for (int j_y = 0; j_y < b.get_height(); j_y++) {
			for (int j_x = 0; j_x < b.get_width(); j_x++) {
				for (int k_y = 0; k_y < filter_weights.get_height(); k_y++) {
					if (k_y + offset_y < j_y - radius_width)
						continue;
					if (k_y + offset_y > j_y + radius_width)
						continue;

					for (int k_x = 0; k_x < filter_weights.get_width(); k_x++) {
						if (k_x + offset_x < j_x - radius_width)
							continue;
						if (k_x + offset_x > j_x + radius_width)
							continue;
						b(j_x, j_y) += filter_weights(k_x, k_y).direct_product(filter_weights(k_x + offset_x - j_x + radius_width, k_y + offset_y - j_y + radius_height));
					}
				}
			}
		}
	}

	vector_fixed<double, 4> b_value(array2d<vector_fixed<double, 4> >& b, const int i_x, const int i_y, const int j_x, const int j_y)
	{
		int radius_width = (b.get_width() - 1) / 2, radius_height = (b.get_height() - 1) / 2;
		int k_x = j_x - i_x + radius_width;
		int k_y = j_y - i_y + radius_height;
		if (k_x < 0 || k_y < 0 || k_x >= b.get_width() || k_y >= b.get_height())
			return vector_fixed<double, 4>();
		return b(k_x, k_y);
	}

	void compute_a_image(const array2d<vector_fixed<double, 4> >& image, array2d<vector_fixed<double, 4> >& b, array2d<vector_fixed<double, 4> >& a)
	{
		int radius_width = (b.get_width() - 1) / 2, radius_height = (b.get_height() - 1) / 2;
		for (int i_y = 0; i_y < a.get_height(); i_y++) {
			for (int i_x = 0; i_x < a.get_width(); i_x++) {
				for (int j_y = i_y - radius_height; j_y <= i_y + radius_height; j_y++) {
					if (j_y < 0)
						j_y = 0;
					if (j_y >= a.get_height())
						break;

					for (int j_x = i_x - radius_width; j_x <= i_x + radius_width; j_x++) {
						if (j_x < 0)
							j_x = 0;
						if (j_x >= a.get_width())
							break;

						a(i_x, i_y) += b_value(b, i_x, i_y, j_x, j_y).direct_product(image(j_x, j_y));
					}
				}
				a(i_x, i_y) *= -2.0;
			}
		}
	}

	void sum_coarsen(const array2d<vector_fixed<double, 4> >& fine, array2d<vector_fixed<double, 4> >& coarse)
	{
		for (int y = 0; y<coarse.get_height(); y++) {
			for (int x = 0; x<coarse.get_width(); x++) {
				coarse(x, y) = fine(x * 2, y * 2);
				if (y * 2 + 1 < fine.get_height())
					coarse(x, y) += fine(x * 2, y * 2 + 1);

				if (x * 2 + 1 < fine.get_width()) {
					coarse(x, y) += fine(x * 2 + 1, y * 2);
					if (y * 2 + 1 < fine.get_height())
						coarse(x, y) += fine(x * 2 + 1, y * 2 + 1);
				}
			}
		}
	}

	template <typename T, int length>
	array2d<T> extract_vector_layer_2d(const array2d<vector_fixed<T, length> >& s, short k)
	{
		array2d<T> result(s.get_width(), s.get_height());
		for (int i = 0; i < s.get_width(); i++) {
			for (int j = 0; j < s.get_height(); j++)
				result(i, j) = s(i, j)[k];
		}
		return result;
	}

	template <typename T, int length>
	vector<T> extract_vector_layer_1d(const vector<vector_fixed<T, length> >& s, short k)
	{
		vector<T> result(s.size());
		for (UINT i = 0; i < s.size(); i++)
			result[i] = s[i][k];

		return result;
	}

	UINT best_match_color(const array3d<double>& vars, const int i_x, const int i_y, const UINT nMaxColor)
	{
		UINT max_v = 0;
		UINT index = i_y * vars.get_width() + i_x;

		auto max_weight = vars(i_x, i_y, max_v);
		for (UINT v = 1; v < nMaxColor; v++) {
			const auto& weight = vars(i_x, i_y, v);
			if (weight > max_weight) {
				max_v = v;
				max_weight = weight;
			}
		}

		return max_v;
	}

	void zoom_double(const array3d<double>& smallVal, array3d<double>& big)
	{
		// Simple scaling of the weights array based on mixing the four
		// pixels falling under each fine pixel, weighted by area.
		// To mix the pixels a little, we assume each fine pixel
		// is 1.2 fine pixels wide and high.
		for (int y = 0; y < big.get_height(); y++) {
			for (int x = 0; x< big.get_width(); x++) {
				double left = max(0.0, (x - 0.1) / 2.0), right = min(smallVal.get_width() - 0.001, (x + 1.1) / 2.0);
				double top = max(0.0, (y - 0.1) / 2.0), bottom = min(smallVal.get_height() - 0.001, (y + 1.1) / 2.0);
				int x_left = (int)floor(left), x_right = (int)floor(right);
				int y_top = (int)floor(top), y_bottom = (int)floor(bottom);
				double area = (right - left) * (bottom - top);
				double top_left_weight = (ceil(left) - left) * (ceil(top) - top) / area;
				double top_right_weight = (right - floor(right)) * (ceil(top) - top) / area;
				double bottom_left_weight = (ceil(left) - left) * (bottom - floor(bottom)) / area;
				double bottom_right_weight = (right - floor(right)) * (bottom - floor(bottom)) / area;
				double top_weight = (right - left) * (ceil(top) - top) / area;
				double bottom_weight = (right - left) * (bottom - floor(bottom)) / area;
				double left_weight = (bottom - top) * (ceil(left) - left) / area;
				double right_weight = (bottom - top) * (right - floor(right)) / area;
				for (int z = 0; z<big.get_depth(); z++) {
					if (x_left == x_right && y_top == y_bottom)
						big(x, y, z) = smallVal(x_left, y_top, z);
					else if (x_left == x_right)
						big(x, y, z) = top_weight * smallVal(x_left, y_top, z) + bottom_weight * smallVal(x_left, y_bottom, z);
					else if (y_top == y_bottom)
						big(x, y, z) = left_weight * smallVal(x_left, y_top, z) + right_weight * smallVal(x_right, y_top, z);
					else {
						big(x, y, z) = top_left_weight * smallVal(x_left, y_top, z) + top_right_weight * smallVal(x_right, y_top, z) +
							bottom_left_weight * smallVal(x_left, y_bottom, z) + bottom_right_weight * smallVal(x_right, y_bottom, z);
					}
				}
			}
		}
	}

	void compute_initial_s(array2d<vector_fixed<double, 4> >& s, const array3d<double>& coarse_variables, array2d<vector_fixed<double, 4> >& b)
	{
		const int length = hasSemiTransparency ? 4 : 3;
		int palette_size = s.get_width();
		int coarse_width = coarse_variables.get_width();
		int coarse_height = coarse_variables.get_height();
		int center_x = (b.get_width() - 1) / 2, center_y = (b.get_height() - 1) / 2;
		auto center_b = b_value(b, 0, 0, 0, 0);
		vector_fixed<double, 4> zero_vector;
		for (int v = 0; v < palette_size; v++) {
			for (int alpha = v; alpha < palette_size; alpha++)
				s(v, alpha) = zero_vector;
		}
		for (int i_y = 0; i_y<coarse_height; i_y++) {
			int max_j_y = min(coarse_height, i_y - center_y + b.get_height());
			for (int i_x = 0; i_x<coarse_width; i_x++) {
				int max_j_x = min(coarse_width, i_x - center_x + b.get_width());
				for (int j_y = max(0, i_y - center_y); j_y<max_j_y; j_y++) {
					for (int j_x = max(0, i_x - center_x); j_x<max_j_x; j_x++) {
						if (i_x == j_x && i_y == j_y)
							continue;
						auto& b_ij = b_value(b, i_x, i_y, j_x, j_y);
						for (int v = 0; v<palette_size; v++) {
							auto v1 = coarse_variables(i_x, i_y, v);
							for (int alpha = v; alpha<palette_size; alpha++) {
								auto mult = v1 * coarse_variables(j_x, j_y, alpha);
								for (byte p = 0; p < length; ++p)
									s(v, alpha)[p] += mult * b_ij[p];
							}
						}
					}
				}
				for (int v = 0; v < palette_size; v++)
					s(v, v) += coarse_variables(i_x, i_y, v) * center_b;
			}
		}
	}

	void update_s(array2d<vector_fixed<double, 4> >& s, const array3d<double>& coarse_variables, array2d<vector_fixed<double, 4> >& b,
		const int j_x, const int j_y, const int alpha, const double delta)
	{
		const int length = hasSemiTransparency ? 4 : 3;
		int palette_size = s.get_width();
		auto coarse_width = coarse_variables.get_width();
		auto coarse_height = coarse_variables.get_height();
		int center_x = (b.get_width() - 1) / 2, center_y = (b.get_height() - 1) / 2;
		int max_i_x = min(coarse_width, j_x + center_x + 1);
		int max_i_y = min(coarse_height, j_y + center_y + 1);
		for (int i_y = max(0, j_y - center_y); i_y < max_i_y; i_y++) {
			for (int i_x = max(0, j_x - center_x); i_x < max_i_x; i_x++) {
				auto delta_b_ij = delta * b_value(b, i_x, i_y, j_x, j_y);
				if (i_x == j_x && i_y == j_y)
					continue;
				for (int v = 0; v <= alpha; v++) {
					auto mult = coarse_variables(i_x, i_y, v);
					for(byte p = 0; p < length; ++p)
						s(v, alpha)[p] += mult * delta_b_ij[p];
				}
				for (int v = alpha; v<palette_size; v++) {
					auto mult = coarse_variables(i_x, i_y, v);
					for (byte p = 0; p < length; ++p)
						s(alpha, v)[p] += mult * delta_b_ij[p];
				}
			}
		}
		s(alpha, alpha) += delta * b_value(b, 0, 0, 0, 0);
	}

	void refine_palette(array2d<vector_fixed<double, 4> >& s, const array3d<double>& coarse_variables,
		const array2d<vector_fixed<double, 4> >& a, vector<vector_fixed<double, 4> >& palette)
	{
		// We only computed the half of S above the diagonal - reflect it
		for (int v = 0; v<s.get_width(); v++) {
			for (int alpha = 0; alpha<v; alpha++)
				s(v, alpha) = s(alpha, v);
		}

		vector<vector_fixed<double, 4> > r(palette.size());
		for (int v = 0; v<palette.size(); v++) {
			for (int i_y = 0; i_y<coarse_variables.get_height(); i_y++) {
				for (int i_x = 0; i_x<coarse_variables.get_width(); i_x++)
					r[v] += coarse_variables(i_x, i_y, v) * a(i_x, i_y);
			}
		}

		const int length = hasSemiTransparency ? 4 : 3;
		for (int k = 0; k < length; k++) {
			auto& S_k = extract_vector_layer_2d(s, k);
			auto& R_k = extract_vector_layer_1d(r, k);
			auto& palette_channel = -1.0 * ((2.0 * S_k).matrix_inverse()) * R_k;
			for (UINT v = 0; v < palette.size(); v++) {
				double val = palette_channel[v];
				if (val < 0.0 || isnan(val))
					val = 0.0;
				else if (val > 1.0)
					val = 1.0;

				palette[v][k] = val;

				if (m_transparentPixelIndex >= 0 && !hasSemiTransparency && k > 1) {
					auto argb = Color::MakeARGB(BYTE_MAX, static_cast<byte>(BYTE_MAX * palette[v][0]), static_cast<byte>(BYTE_MAX * palette[v][1]), static_cast<byte>(BYTE_MAX * palette[v][2]));
					if (Color(argb).ToCOLORREF() == Color(m_transparentColor).ToCOLORREF())
						swap(palette[0], palette[v]);
				}
			}			
		}
	}

	void compute_initial_j_palette_sum(array2d<vector_fixed<double, 4> >& j_palette_sum, const array3d<double>& coarse_variables, const vector<vector_fixed<double, 4> >& palette)
	{
		for (int j_y = 0; j_y<coarse_variables.get_height(); j_y++) {
			for (int j_x = 0; j_x<coarse_variables.get_width(); j_x++) {
				vector_fixed<double, 4> palette_sum;
				for (UINT alpha = 0; alpha < palette.size(); alpha++)
					palette_sum += coarse_variables(j_x, j_y, alpha) * palette[alpha];
				j_palette_sum(j_x, j_y) = palette_sum;
			}
		}
	}

	bool spatial_color_quant(const array2d<vector_fixed<double, 4> >& image, array2d<vector_fixed<double, 4> >& filter_weights,
		array2d<UINT>& quantized_image, vector<vector_fixed<double, 4> >& palette,
		const double initial_temperature = 1.0, const double final_temperature = 0.001, const int temps_per_level = 3, const int repeats_per_temp = 1)
	{
		const int length = hasSemiTransparency ? 4 : 3;
		UINT nMaxColor = palette.size();
		int max_coarse_level = //1;
			compute_max_coarse_level(image.get_width(), image.get_height());
		auto p_coarse_variables = new array3d<double>(
			image.get_width() >> max_coarse_level,
			image.get_height() >> max_coarse_level,
			palette.size());
		auto pp_coarse_variables = unique_ptr<array3d<double> >(p_coarse_variables);
		// For syntactic convenience
		auto& coarse_variables = *p_coarse_variables;
		coarse_variables.fill_random();

		double temperature = initial_temperature;

		// Compute a_i, b_{ij} according to (11)
		int extended_neighborhood_width = filter_weights.get_width() * 2 - 1;
		int extended_neighborhood_height = filter_weights.get_height() * 2 - 1;
		array2d<vector_fixed<double, 4> > b0(extended_neighborhood_width,
			extended_neighborhood_height);
		compute_b_array(filter_weights, b0);

		array2d<vector_fixed<double, 4> > a0(image.get_width(), image.get_height());
		compute_a_image(image, b0, a0);

		// Compute a_I^l, b_{IJ}^l according to (18)
		vector<array2d<vector_fixed<double, 4> > > a_vec, b_vec;
		a_vec.emplace_back(a0);
		b_vec.emplace_back(b0);

		int radius_width = (filter_weights.get_width() - 1) / 2, radius_height = (filter_weights.get_height() - 1) / 2;
		int coarse_level;
		for (coarse_level = 1; coarse_level <= max_coarse_level; coarse_level++) {
			array2d<vector_fixed<double, 4> > bi(max(length, b_vec.back().get_width() - 2), max(length, b_vec.back().get_height() - 2));

			for (int J_y = 0; J_y<bi.get_height(); J_y++) {
				for (int J_x = 0; J_x<bi.get_width(); J_x++) {
					for (int i_y = radius_height * 2; i_y<radius_height * 2 + 2; i_y++) {
						for (int i_x = radius_width * 2; i_x<radius_width * 2 + 2; i_x++) {
							for (int j_y = J_y * 2; j_y<J_y * 2 + 2; j_y++) {
								for (int j_x = J_x * 2; j_x<J_x * 2 + 2; j_x++)
									bi(J_x, J_y) += b_value(b_vec.back(), i_x, i_y, j_x, j_y);
							}
						}
					}
				}
			}
			b_vec.emplace_back(bi);

			array2d<vector_fixed<double, 4> > ai(image.get_width() >> coarse_level, image.get_height() >> coarse_level);
			sum_coarsen(a_vec.back(), ai);
			a_vec.emplace_back(ai);
		}

		// Multiscale annealing
		coarse_level = max_coarse_level;
		const int iters_per_level = temps_per_level;
		double temperature_multiplier = pow(final_temperature / initial_temperature, 1.0 / (max(length, max_coarse_level * iters_per_level)));

		int iters_at_current_level = 0;
		bool skip_palette_maintenance = false;
		array2d<vector_fixed<double, 4> > s(palette.size(), palette.size());
		compute_initial_s(s, coarse_variables, b_vec[coarse_level]);
		auto p_palette_sum = make_unique<array2d< vector_fixed<double, 4> > >(coarse_variables.get_width(), coarse_variables.get_height());
		auto j_palette_sum = p_palette_sum.get();
		compute_initial_j_palette_sum(*j_palette_sum, coarse_variables, palette);

		while (coarse_level >= 0 || temperature > final_temperature) {
			// Need to reseat this reference in case we changed p_coarse_variables
			auto& coarse_variables = *p_coarse_variables;
			auto& a = a_vec[coarse_level];
			auto& b = b_vec[coarse_level];
			vector_fixed<double, 4> middle_b = b_value(b, 0, 0, 0, 0);

			int center_x = (b.get_width() - 1) / 2, center_y = (b.get_height() - 1) / 2;
			int step_counter = 0;
			for (int repeat = 0; repeat < repeats_per_temp; repeat++) {
				int pixels_changed = 0, pixels_visited = 0;
				deque<pair<int, int> > visit_queue;
				random_permutation_2d(coarse_variables.get_width(), coarse_variables.get_height(), visit_queue);

				// Compute 2*sum(j in extended neighborhood of i, j != i) b_ij

				while (!visit_queue.empty()) {
					// If we get to 10% above initial size, just revisit them all
					if ((int)visit_queue.size() > coarse_variables.get_width() * coarse_variables.get_height() * 1.1) {
						visit_queue.clear();
						random_permutation_2d(coarse_variables.get_width(), coarse_variables.get_height(), visit_queue);
					}

					int i_x = visit_queue.front().first, i_y = visit_queue.front().second;
					visit_queue.pop_front();

					// Compute (25)
					vector_fixed<double, 4> p_i;
					for (int y = 0; y<b.get_height(); y++) {
						int j_y = y - center_y + i_y;
						if (j_y < 0 || j_y >= coarse_variables.get_height())
							continue;
						for (int x = 0; x<b.get_width(); x++) {
							int j_x = x - center_x + i_x;
							if (i_x == j_x && i_y == j_y)
								continue;
							if (j_x < 0 || j_x >= coarse_variables.get_width())
								continue;
							auto& b_ij = b_value(b, i_x, i_y, j_x, j_y);
							auto& j_pal = (*j_palette_sum)(j_x, j_y);
							for (byte p = 0; p < length; ++p)
								p_i[p] += b_ij[p] * j_pal[p];
						}
					}
					p_i *= 2.0;
					p_i += a(i_x, i_y);

					vector<double> meanfield_logs, meanfields;
					double max_meanfield_log = -numeric_limits<double>::infinity();
					double meanfield_sum = 0.0;
					for (UINT v = 0; v < palette.size(); v++) {
						// Update m_{pi(i)v}^I according to (23)
						// We can subtract an arbitrary factor to prevent overflow,
						// since only the weight relative to the sum matters, so we
						// will choose a value that makes the maximum e^100.
						meanfield_logs.push_back(-(palette[v].dot_product(p_i + middle_b.direct_product(palette[v]))) / temperature);
						if (meanfield_logs.back() > max_meanfield_log)
							max_meanfield_log = meanfield_logs.back();
					}
					for (UINT v = 0; v < palette.size(); v++) {
						meanfields.push_back(exp(meanfield_logs[v] - max_meanfield_log + 100));
						meanfield_sum += meanfields.back();
					}
					if (meanfield_sum == 0)
						return false;

					auto old_max_v = best_match_color(coarse_variables, i_x, i_y, palette.size());
					auto& j_pal = (*j_palette_sum)(i_x, i_y);
					for (UINT v = 0; v < palette.size(); v++) {
						double new_val = meanfields[v] / meanfield_sum;
						// Prevent the matrix S from becoming singular
						if (new_val <= 0)
							new_val = 1e-10;
						if (new_val >= 1)
							new_val = 1 - 1e-10;
						double delta_m_iv = new_val - coarse_variables(i_x, i_y, v);
						UINT index = i_y * coarse_variables.get_width() + i_x;

						coarse_variables(i_x, i_y, v) = new_val;
						for (byte p = 0; p < length; ++p)
							j_pal[p] += delta_m_iv * palette[v][p];

						if (abs(delta_m_iv) > 0.001 && !skip_palette_maintenance)
							update_s(s, coarse_variables, b, i_x, i_y, v, delta_m_iv);
					}
					auto max_v = best_match_color(coarse_variables, i_x, i_y, palette.size());
					// Only consider it a change if the colors are different enough
					if ((palette[max_v] - palette[old_max_v]).norm_squared() >= 1.0 / (255.0 * 255.0)) {
						pixels_changed++;
						// We don't add the outer layer of pixels , because
						// there isn't much weight there, and if it does need
						// to be visited, it'll probably be added when we visit
						// neighboring pixels.
						// The commented out loops are faster but cause a little bit of distortion
						//for (int y=center_y-1; y<center_y+1; y++) {
						//   for (int x=center_x-1; x<center_x+1; x++) {
						for (int y = min(1, center_y - 1); y < max(b.get_height() - 1, center_y + 1); y++) {
							int j_y = y - center_y + i_y;
							if (j_y < 0 || j_y >= coarse_variables.get_height())
								continue;
							for (int x = min(1, center_x - 1); x < max(b.get_width() - 1, center_x + 1); x++) {
								int j_x = x - center_x + i_x;
								if (j_x < 0 || j_x >= coarse_variables.get_width())
									continue;
								visit_queue.emplace_back(j_x, j_y);
							}
						}
					}
					pixels_visited++;

					// Show progress with dots - in a graphical interface,
					// we'd show progressive refinements of the image instead,
					// and maybe a palette preview.
					step_counter++;
				}
				if (skip_palette_maintenance)
					compute_initial_s(s, *p_coarse_variables, b_vec[coarse_level]);

				refine_palette(s, coarse_variables, a, palette);
				compute_initial_j_palette_sum(*j_palette_sum, coarse_variables, palette);
			}

			iters_at_current_level++;
			skip_palette_maintenance = false;
			if ((temperature <= final_temperature || coarse_level > 0) && iters_at_current_level >= iters_per_level) {
				if (--coarse_level < 0)
					break;
				auto p_new_coarse_variables = new array3d<double>(image.get_width() >> coarse_level, image.get_height() >> coarse_level, palette.size());
				zoom_double(coarse_variables, *p_new_coarse_variables);
				pp_coarse_variables.reset(p_new_coarse_variables);
				p_coarse_variables = pp_coarse_variables.get();
				iters_at_current_level = 0;
				p_palette_sum.reset(new array2d<vector_fixed<double, 4> >(p_coarse_variables->get_width(), p_coarse_variables->get_height()));
				j_palette_sum = p_palette_sum.get();
				compute_initial_j_palette_sum(*j_palette_sum, *p_coarse_variables, palette);
				skip_palette_maintenance = true;
			}
			if (temperature > final_temperature)
				temperature *= temperature_multiplier;
		}

		for (int i_x = 0; i_x < image.get_width(); i_x++) {
			for (int i_y = 0; i_y < image.get_height(); i_y++)
				quantized_image(i_x, i_y) = best_match_color(*p_coarse_variables, i_x, i_y, palette.size());
		}

		return true;
	}

	bool ProcessImagePixels(Bitmap* pDest, const ColorPalette* pPalette, array2d<UINT>& quantized_image)
	{
		pDest->SetPalette(pPalette);

		BitmapData targetData;
		UINT w = pDest->GetWidth();
		UINT h = pDest->GetHeight();

		Status status = pDest->LockBits(&Gdiplus::Rect(0, 0, w, h), ImageLockModeWrite, pDest->GetPixelFormat(), &targetData);
		if (status != Ok) {
			cerr << "Cannot write image" << endl;
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

	bool SpatialQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		if (nMaxColors > 256)
			nMaxColors = 256;
		
		const UINT bitDepth = GetPixelFormatSize(pSource->GetPixelFormat());
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		hasSemiTransparency = false;
		m_transparentPixelIndex = -1;
		int pixelIndex = 0;
		array2d<vector_fixed<double, 4> > pixels(bitmapWidth, bitmapHeight);
		if (bitDepth <= 16) {
			for (UINT y = 0; y < bitmapHeight; y++) {
				for (UINT x = 0; x < bitmapWidth; x++, pixelIndex++) {
					Color color;
					pSource->GetPixel(x, y, &color);
					pixels(x, y)[0] = color.GetR() / ((double)BYTE_MAX);
					pixels(x, y)[1] = color.GetG() / ((double)BYTE_MAX);
					pixels(x, y)[2] = color.GetB() / ((double)BYTE_MAX);
					pixels(x, y)[3] = color.GetA() / ((double)BYTE_MAX);

					if (color.GetA() < BYTE_MAX) {
						hasSemiTransparency = true;
						m_transparentPixelIndex = pixelIndex;
						if (color.GetA() == 0)
							m_transparentColor = color.GetValue();
					}

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

				for (UINT x = 0; x < bitmapWidth; x++, pixelIndex++) {	// ...for each pixel...
					byte pixelBlue = *pPixelSource++;
					byte pixelGreen = *pPixelSource++;
					byte pixelRed = *pPixelSource++;
					byte pixelAlpha = bitDepth < 32 ? BYTE_MAX : *pPixelSource++;

					pixels(x, y)[0] = pixelRed / ((double)BYTE_MAX);
					pixels(x, y)[1] = pixelGreen / ((double)BYTE_MAX);
					pixels(x, y)[2] = pixelBlue / ((double)BYTE_MAX);
					pixels(x, y)[3] = pixelAlpha / ((double)BYTE_MAX);

					if (pixelAlpha < BYTE_MAX) {
						hasSemiTransparency = true;
						m_transparentPixelIndex = pixelIndex;
						if (pixelAlpha == 0)
							m_transparentColor = Color::MakeARGB(pixelAlpha, pixelRed, pixelGreen, pixelBlue);
					}
				}

				pRowSource += strideSource;
			}

			pSource->UnlockBits(&data);
		}

		const int length = hasSemiTransparency ? 4 : 3;
		double dithering_level = 1.0;
		array2d<vector_fixed<double, 4> > filter3_weights(3, 3);
		double sum = 0.0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double value = exp(-sqrt((double)((i - 2)*(i - 2) + (j - 2)*(j - 2))) / (dithering_level * dithering_level));
				for (short k = 0; k < length; k++)
					sum += filter3_weights(i, j)[k] = value;
			}
		}
		sum /= 3.0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (short k = 0; k < length; k++)
					filter3_weights(i, j)[k] /= sum;
			}
		}

		array2d<UINT> quantized_image(bitmapWidth, bitmapHeight);
		vector<vector_fixed<double, 4> > palette(nMaxColors);
		for (UINT i = 0; i<nMaxColors; i++) {
			for (byte p = 0; p < length; ++p)
				palette[i][p] = ((double)rand()) / RAND_MAX;
		}

		if (!spatial_color_quant(pixels, filter3_weights, quantized_image, palette))
			return false;

		auto pPaletteBytes = make_unique<byte[]>(pDest->GetPaletteSize());
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors > 2) {
			/* Fill palette */
			for (UINT k = 0; k<nMaxColors; ++k)
				pPalette->Entries[k] = Color::MakeARGB(hasSemiTransparency ? static_cast<byte>(BYTE_MAX * palette[k][3]) : BYTE_MAX, static_cast<byte>(BYTE_MAX * palette[k][0]), static_cast<byte>(BYTE_MAX * palette[k][1]), static_cast<byte>(BYTE_MAX * palette[k][2]));
			
			if (m_transparentPixelIndex >= 0) {
				UINT k = quantized_image[m_transparentPixelIndex];
				if (nMaxColors > 2)
					pPalette->Entries[k] = m_transparentColor;
				else if (pPalette->Entries[k] != m_transparentColor)
					swap(pPalette->Entries[0], pPalette->Entries[1]);
			}
		}
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette->Entries[0] = m_transparentColor;
				pPalette->Entries[1] = Color::Black;
			}
			else {
				pPalette->Entries[1] = Color::Black;
				pPalette->Entries[0] = Color::White;
			}
		}

		return ProcessImagePixels(pDest, pPalette, quantized_image);
	}

}