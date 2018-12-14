#pragma once
#include <memory>
#include <vector>

using namespace std;

namespace EdgeAwareSQuant
{
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

		double dot_product(const vector_fixed<T, length>& rhs) {
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
		array2d()
		{
			this->width = 0;
			this->height = 0;
		}

		array2d(int width, int height)
		{
			reset(width, height);
		}

		array2d(const array2d<T>& rhs)
		{
			width = rhs.get_width();
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

		void reset(int width, int height)
		{
			this->width = width;
			this->height = height;
			data = make_unique<T[]>(width * height);
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
			if (mult != 0) {
				for (int i = 0; i<get_width(); i++)
					(*this)(i, to_row) += mult * (*this)(i, from_row);
			}

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
				auto detA = a(i, i);
				float val = (detA != 0.0f) ? 1.0f / detA : 0.0f;
				result.multiply_row_scalar(i, val);
				multiply_row_scalar(i, val);
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
	class Mat
	{
	public:
		Mat()
		{
			this->width = 0;
			this->height = 0;
		}

		Mat(int height, int width)
		{
			reset(height, width);
		}

		Mat(const Mat<T>& rhs)
		{
			width = rhs.get_width();
			height = rhs.height;
			data = make_unique<T[]>(width * height);
			for (int i = 0; i<(width * height); i++)
				data[i] = rhs.data[i];
		}

		inline T& at(int row, int col)
		{
			return data[row * width + col];
		}
		
		inline const T& at(int row, int col) const
		{
			return data[row * width + col];
		}

		inline T& operator()(int row, int col)
		{
			return data[row * width + col];
		}
		
		inline const T& operator()(int row, int col) const
		{
			return data[row * width + col];
		}

		inline int get_width() const { return width; }
		inline int get_height() const { return height; }

		Mat<T>& operator*=(const T scalar) {
			for (int i = 0; i<(width * height); i++)
				data[i] *= scalar;
			return *this;
		}

		Mat<T>& operator/=(const T scalar) {
			for (int i = 0; i<(width * height); i++)
				data[i] /= scalar;
			return *this;
		}

		Mat<T> operator*(const T scalar) {
			Mat<T> result(*this);
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

		void reset(int height, int width)
		{
			this->width = width;
			this->height = height;
			data = make_unique<T[]>(width * height);
		}

	private:
		unique_ptr<T[]> data;
		int width, height;
	};

	template <typename T>
	Mat<T> operator*(T scalar, Mat<T> a) {
		return a * scalar;
	}

	// =============================================================
	// Quantizer objects and functions
	//
	// COVERED CODE IS PROVIDED UNDER THIS LICENSE ON AN "AS IS" BASIS, WITHOUT WARRANTY
	// OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
	// THAT THE COVERED CODE IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE
	// OR NON-INFRINGING. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE COVERED
	// CODE IS WITH YOU. SHOULD ANY COVERED CODE PROVE DEFECTIVE IN ANY RESPECT, YOU (NOT
	// THE INITIAL DEVELOPER OR ANY OTHER CONTRIBUTOR) ASSUME THE COST OF ANY NECESSARY
	// SERVICING, REPAIR OR CORRECTION. THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL
	// PART OF THIS LICENSE. NO USE OF ANY COVERED CODE IS AUTHORIZED HEREUNDER EXCEPT UNDER
	// THIS DISCLAIMER.
	//
	// Use at your own risk!
	// =============================================================

	class EdgeAwareSQuantizer
	{
		public:
			bool QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither = true);
	};
}