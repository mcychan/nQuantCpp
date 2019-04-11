#pragma once
/*
 
 Copyright (c) 2015, M. Emre Celebi
 Copyright (c) 2019 Miller Cy Chan
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */

#include "stdafx.h"
#include "DivQuantizer.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include <algorithm>
#include <unordered_map>
#include <type_traits>

namespace DivQuant
{
	const int COLOR_HASH_SIZE = 20023;
	bool hasSemiTransparency = false;
	int m_transparentPixelIndex = -1;
	ARGB m_transparentColor = Color::Transparent;
	unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;

	struct Bucket
	{
		byte value = 0;
		ARGB argb = Color::Transparent;
		shared_ptr<Bucket> next;
	};
	
	template <
		typename T, //real type
		typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type
	> struct Pixel
	{
		T alpha = 0;
		double L, A, B;
		ARGB argb;
		T weight = 0;
	};

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

	// This method will dedup unique pixels and subsample pixels
	// based on dec_factor. When dec_factor is 1 then this method
	// would not do anything if the input is already unique, use
	// unique_colors_as_doubles() in that case.
	unique_ptr<double[]> calc_color_table(const ARGB* inPixels,
                  const UINT numPixels,
                  ARGB* outPixels,
                  const UINT numRows,
                  const UINT numCols,
                  const int dec_factor,
                  UINT& num_colors)
	{ 
		unique_ptr<double[]> weights;
		if(dec_factor <= 0) {
			cerr << "Decimation factor ( " << dec_factor << " ) should be positive !\n";    
			return weights;
		}
  
		vector<shared_ptr<Bucket> > hash_table(COLOR_HASH_SIZE);
		num_colors = 0;
  
		for (UINT ir = 0; ir < numRows; ir += dec_factor) {
			for (UINT ic = 0; ic < numCols; ic += dec_factor) {
				Color c(inPixels[ic + (ir * numRows)]);				
      
				/* Determine the bucket */
				int hash = getARGBIndex(c, hasSemiTransparency, m_transparentPixelIndex) % COLOR_HASH_SIZE;				
				shared_ptr<Bucket> bucket;
				/* Search for the color in the bucket chain */
				for (bucket = hash_table[hash]; bucket.get() != nullptr; bucket = bucket->next) {
					if (bucket->argb == c.GetValue()) {
						/* This color exists in the hash table */
						break;
					}
				}
	  
				if (bucket.get() != nullptr) {
					/* This color exists in the hash table */
					bucket->value++;
				}
				else {
					++num_colors;

					/* Create a new bucket entry for this color */
					bucket = make_shared<Bucket>();					
					bucket->value = 1;
					bucket->argb = c.GetValue();
					bucket->next = hash_table[hash];
					hash_table[hash] = bucket;
				}
			}
		}
  
		weights = make_unique<double[]>(num_colors);
  
		/* Normalization factor to obtain color frequencies to color probabilities */
		double norm_factor = 1.0 / (ceil(numRows / (double) dec_factor) * ceil(numCols / (double) dec_factor));
  
		for (int index = 0, hash = 0; hash < COLOR_HASH_SIZE; ++hash) {
			for (auto bucket = hash_table[hash]; bucket.get() != nullptr; bucket = bucket->next) {
				outPixels[index] = bucket->argb;      
				weights[index++] = norm_factor * bucket->value;
			}
		}
  
		hash_table.clear();
		return weights;
	}

	double get_double_scale(const UINT numPixels)
	{
		const int numRows = 1;
		const int dec_factor = 1;

		return 1.0 / (ceil(numRows / (double) dec_factor) * ceil(numPixels / (double) dec_factor));
	}

	static inline bool asc_weighted_pixel(const Pixel<int>& a, const Pixel<int>& b) {
		return a.weight < b.weight;
	}

	static void sort_color(Pixel<int>* cmap, const int num_colors)
	{
		vector<Pixel<int> > pixelVec(num_colors);
		for (int i = 0; i < num_colors; ++i)
			pixelVec[i] = cmap[i];

		sort(begin(pixelVec), end(pixelVec), asc_weighted_pixel);
		for (int i = 0; i < num_colors; ++i)
			cmap[i] = pixelVec[i];
	}

	void map_colors_mps(const ARGB* inPixelsPtr, UINT numPixels, short* qPixels, ColorPalette* pPalette)
	{
		const UINT colormapSize = pPalette->Count;
		const int size_lut_init = 4 * BYTE_MAX + 1;
		const int max_sum = 4 * BYTE_MAX;
  
		auto lut_init = make_unique<int[]>(size_lut_init);
	  
		auto cmap = make_unique<Pixel<int>[]>(colormapSize);
		for (UINT i = 0; i < colormapSize; ++i) {
			Color c(pPalette->Entries[i]);
			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
			
			auto& pi = cmap[i];
			pi.alpha = c.GetA();
			pi.L = lab1.L;
			pi.A = lab1.A;
			pi.B = lab1.B;
			pi.argb = c.GetValue();
			pi.weight = c.GetA() + c.GetR() + c.GetG() + c.GetB();
		}
	  
		const int size_lut_ssd = 2 * max_sum + 1;
		auto lut_ssd_buffer = make_unique<int[]>(size_lut_ssd);
	  
		auto lut_ssd = lut_ssd_buffer.get() + max_sum;
		lut_ssd[0] = 0;
	  
		// Premultiply the LUT entries by (1/4) -- see below
		for (int ik = 1; ik <= max_sum; ++ik)
			lut_ssd[-ik] = lut_ssd[ik] = (int) (sqr(ik) / 4.0);
	  
		// Sort the palette by the sum of color components.	  
		sort_color(cmap.get(), colormapSize);

		for (UINT i = 0; i < colormapSize; ++i)
			pPalette->Entries[i] = cmap[i].argb;
	  
		// Calculate the LUT
		int low = (colormapSize >= 2) ? (int) (0.5 * (cmap[0].weight + cmap[1].weight) + 0.5) : 1;
	  
		int high = (colormapSize >= 2) ? (int) (0.5 * (cmap[colormapSize - 2].weight + cmap[colormapSize - 1].weight) + 0.5) : 1;
	  
		for (int ik = high; ik < size_lut_init; ++ik)
			lut_init[ik] = colormapSize - 1;
	  
		for (UINT ic = 1; ic < colormapSize - 1; ++ic) {
			low = (int) (0.5 * (cmap[ic - 1].weight + cmap[ic].weight) + 0.5);  // round
			high = (int) (0.5 * (cmap[ic].weight + cmap[ic + 1].weight) + 0.5); // round
		
			for (int ik = low; ik < high; ++ik)
				lut_init[ik] = ic;
		}
	  
		for (UINT ik = 0; ik < numPixels; ++ik) {
			Color c(inPixelsPtr[ik]);
			int sum = c.GetA() + c.GetR() + c.GetG() + c.GetB();
		
			// Determine the initial searched colour cinit in the palette for cp.
			int index = lut_init[sum];		
		
			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
			// Calculate the squared Euclidean distance between cp and cinit
			UINT min_dist = abs(c.GetA() - cmap[index].alpha) + abs(lab1.L - cmap[index].L) + abs(lab1.A - cmap[index].A) + abs(lab1.B - cmap[index].B);
			int upi = index, downi = index;
			bool up = true, down = true;
			while (up || down) {
				if (up) {				
					if (++upi > (colormapSize - 1) || lut_ssd[sum - cmap[upi].weight] >= min_dist) {
						// Terminate the search in UP direction
						up = false;          
					}
					else {
						UINT dist = abs(c.GetA() - cmap[upi].alpha) + abs(lab1.L - cmap[upi].L) + abs(lab1.A - cmap[upi].A) + abs(lab1.B - cmap[upi].B);          
						if (dist < min_dist) {
							min_dist = dist;
							index = upi;            
						}
					}
				}
			  
				if (down) {				
					if (--downi < 0 || lut_ssd[sum - cmap[downi].weight] >= min_dist) {
						// Terminate the search in DOWN direction
						down = false;
					}
					else {
						UINT dist = abs(c.GetA() - cmap[downi].alpha) + abs(lab1.L - cmap[downi].L) + abs(lab1.A - cmap[downi].A) + abs(lab1.B - cmap[downi].B);
						if (dist < min_dist) {
							min_dist = dist;
							index = downi;            
						}
					}
				}
			}

			qPixels[ik] = index;
		}
	}

	// MT  : type of the member attribute, either byte or UINT
	template <typename MT>
	void DivQuantClusterInitMeanAndVar(const int num_points, const ARGB* data, const double data_weight, double* weightsPtr, Pixel<double>& total_mean, Pixel<double>& total_var)
	{
		double mean_alpha = 0.0, mean_L = 0.0, mean_A = 0.0, mean_B = 0.0;
		double var_alpha = 0.0, var_L = 0.0, var_A = 0.0, var_B = 0.0;
  
		for (int ip = 0; ip < num_points; ++ip) {
			Color c(data[ip]);
			CIELABConvertor::Lab lab1;
			getLab(c, lab1);
    
			if (weightsPtr == nullptr) {
				mean_alpha += c.GetA();
				mean_L += lab1.L;
				mean_A += lab1.A;
				mean_B += lab1.B;
      
				var_alpha += sqr(c.GetA());
				var_L += sqr(lab1.L);
				var_A += sqr(lab1.A);
				var_B += sqr(lab1.B);
			} else {
				// non-uniform weights      
				double tmp_weight = weightsPtr[ip];
      
				mean_alpha += tmp_weight * c.GetA();
				mean_L += tmp_weight * lab1.L;
				mean_A += tmp_weight * lab1.A;
				mean_B += tmp_weight * lab1.B;
      
				var_alpha += tmp_weight * sqr(c.GetA());
				var_L += tmp_weight * sqr(lab1.L);
				var_A += tmp_weight * sqr(lab1.A);
				var_B += tmp_weight * sqr(lab1.B);
			}
		}
  
		if (weightsPtr == nullptr) {
			// In uniform weight case do the multiply outside the loop
			mean_alpha *= data_weight;
			mean_L *= data_weight;
			mean_A *= data_weight;
			mean_B *= data_weight;
    
			var_alpha *= data_weight;
			var_L *= data_weight;
			var_A *= data_weight;
			var_B *= data_weight;
		}
  
		var_alpha -= sqr(mean_alpha);
		var_L -= sqr(mean_L);
		var_A -= sqr(mean_A);
		var_B -= sqr(mean_B);
  
		// Copy data to user supplied pointers
		total_mean.alpha = mean_alpha;
		total_mean.L = mean_L;
		total_mean.A = mean_A;
		total_mean.B = mean_B;
  
		total_var.alpha = var_alpha;
		total_var.L = var_L;
		total_var.A = var_A;
		total_var.B = var_B;
	}

	// This method defines a clustering approach that divides the input into
	// roughly equally sized clusters until N clusters is reached or the
	// clusters can be divided no more.

	// MT  : type of the member attribute, either byte or UINT
	template <typename MT>
	void DivQuantCluster(const int num_points, ARGB* data, ARGB* tmp_buffer, const double data_weight, double* weightsPtr,
        const int num_bits, const int max_iters, ColorPalette* pPalette, UINT& nMaxColors)
	{  
		const UINT num_colors = nMaxColors;
	  
		bool apply_lkm = 0 < max_iters; /* indicates whether or not LKM is to be applied */
		int max_iters_m1 = max_iters - 1;  
		
		int tmp_buffer_used = 0; // Capacity in num points that can be stored in tmp_data

		// The member array is either byte or UINT.
		auto member = make_unique<MT[]>(num_points);

	  
		unique_ptr<int[]> point_index;
	  
		auto weight = make_unique<double[]>(num_colors); /* total weight of each cluster */
	  
		/*
		* Contains the size of each cluster. The size of a cluster is
		* actually the number unique colors that it represents.
		*/
		auto size = make_unique<int[]>(num_colors);
	  
		auto tse = make_unique<double[]>(num_colors); /* total squared error of each cluster */
	  
		auto mean = make_unique<Pixel<double>[]>(num_colors); /* componentwise mean (centroid) of each cluster */
	  
		auto var = make_unique<Pixel<double>[]>(num_colors); /* componentwise variance of each cluster */		
	  
		/* Cluster 0 is always the first cluster to be split */
		int old_index = 0; /* index of C or C1 */
	  
		/* First cluster to be split contains the entire data set */
		weight[old_index] = 1.0;
	  
		int tmp_num_points = num_points; /* number of points in C */
		/*
		* # points is not the same as # pixels. Each point represents
		* potentially multiple pixels with a specific color.
		*/	
		size[old_index] = tmp_num_points;
	  
		/* Perform ( NUM_COLORS - 1 ) splits */
		/*
		OLD_INDEX denotes the index of the cluster to be split.
		When cluster OLD_INDEX is split, the indexes of the two subclusters
		are given by OLD_INDEX and NEW_INDEX, respectively.
		*/
		UINT new_index = 1; /* index of C2 */
		int new_size = 0; /* size of C2 */
		double cut_pos; /* cutting position */
		double total_weight; /* weight of C */
		double max_val;
		
		Pixel<double> total_mean; // componentwise mean of C
		Pixel<double> total_var; // componentwise variance of C
  
		double new_weight; /* weight of C2 */
		double lhs;
		auto tmp_data = data; /* temporary data set (holds the cluster to be split) */
		double tmp_weight; /* weight of a particular pixel */ 
		for (; new_index < num_colors; ++new_index) {
			/* STEPS 1 & 2: DETERMINE THE CUTTING AXIS AND POSITION */    
			total_weight = weight[old_index];
		
			if (new_index == 1)
				DivQuantClusterInitMeanAndVar<MT>(num_points, data, data_weight, weightsPtr, total_mean, total_var);
			else {
				// Cluster mean/variance has already been calculated
				total_mean.alpha = mean[old_index].alpha;
				total_mean.L = mean[old_index].L;
				total_mean.A = mean[old_index].A;
				total_mean.B = mean[old_index].B;
			  
				total_var.alpha = var[old_index].alpha;
				total_var.L = var[old_index].L;
				total_var.A = var[old_index].A;
				total_var.B = var[old_index].B;    
			}
		
			/* Determine the axis with the greatest variance */    
			max_val = total_var.alpha;
			byte cut_axis = 0; /* index of the cutting axis */
			cut_pos = total_mean.alpha;
			if (max_val < total_var.L) {
				max_val = total_var.L;
				cut_axis = 1;
				cut_pos = total_mean.L;
			}
			
			if (max_val < total_var.A) {
				max_val = total_var.A;
				cut_axis = 2;
				cut_pos = total_mean.A;
			}
		
			if (max_val < total_var.B) {
				max_val = total_var.B;
				cut_axis = 3;
				cut_pos = total_mean.B;
			}
		
			auto& new_mean = mean[new_index]; /* componentwise mean of C2 */
			auto& new_var = var[new_index]; /* componentwise variance of C2 */
		
			// Reset the statistics of the new cluster
			new_weight = 0.0;
			UINT new_weight_count = 0;
			new_mean.alpha = new_mean.L = new_mean.A = new_mean.B = 0.0;
		
			if (!apply_lkm) {
				new_size = 0;
				new_var.alpha = new_var.L = new_var.A = new_var.B = 0.0;
			}
			
			// STEP 3: SPLIT THE CLUSTER OLD_INDEX    
			for (int ip = 0; ip < tmp_num_points; ) {
				double new_mean_alpha = 0, new_mean_L = 0, new_mean_A = 0, new_mean_B = 0;
		  
				double new_var_alpha = 0, new_var_L = 0, new_var_A = 0, new_var_B = 0;
		  
				int maxLoopOffset = 0xFFFF;
				int numLeft = (tmp_num_points - ip);
				if (numLeft < maxLoopOffset)
					maxLoopOffset = numLeft;

				maxLoopOffset += ip;
		  
				for (; ip < maxLoopOffset; ++ip) {
					Color c(tmp_data[ip]);
					CIELABConvertor::Lab lab1;
					getLab(c, lab1);
			
					double proj_val = c.GetA();
					if(cut_axis == 1)
						proj_val = lab1.L;
					else if(cut_axis == 2)
						proj_val = lab1.A;
					else if(cut_axis == 3)
						proj_val = lab1.B; /* projection of a data point on the cutting axis */
			
					if (cut_pos < proj_val) {          
						if (weightsPtr == nullptr) {
							new_mean_alpha += c.GetA();
							new_mean_L += lab1.L;
							new_mean_A += lab1.A;
							new_mean_B += lab1.B;
						}
						else {
							// non-uniform weights            
							int pointindex = ip;
							if (point_index.get())
								pointindex = point_index[ip];

							tmp_weight = weightsPtr[pointindex];
				
							new_mean.alpha += tmp_weight * c.GetA();
							new_mean.L += tmp_weight * lab1.L;
							new_mean.A += tmp_weight * lab1.A;
							new_mean.B += tmp_weight * lab1.B;
						}
			  
						// Update the point membership and variance/size of the new cluster
						if (!apply_lkm) {
							int pointindex = ip;
							if (point_index.get())
								pointindex = point_index[ip];

							member[pointindex] = new_index;
				
							if (weightsPtr == nullptr) {
								new_var_alpha += sqr(c.GetA());
								new_var_L += sqr(lab1.L);
								new_var_A += sqr(lab1.A);
								new_var_B += sqr(lab1.B);
							}
							else {
								// non-uniform weights              
								// tmp_weight already set above in loop
								new_var.alpha += tmp_weight * sqr(c.GetA());
								new_var.L += tmp_weight * sqr(lab1.L);
								new_var.A += tmp_weight * sqr(lab1.A);
								new_var.B += tmp_weight * sqr(lab1.B);
							}
				
							++new_size;
						}
			  
						// Update the weight of the new cluster          
						if (weightsPtr == nullptr)
							++new_weight_count;
						else
							new_weight += tmp_weight;
					}        
				} // end foreach tmp_num_points inner loop
		  
				if (weightsPtr == nullptr) {
					new_mean.alpha += new_mean_alpha;
					new_mean.L += new_mean_L;
					new_mean.A += new_mean_A;
					new_mean.B += new_mean_B;
			
					if (!apply_lkm) {
						new_var.alpha += new_var_alpha;
						new_var.L += new_var_L;
						new_var.A += new_var_A;
						new_var.B += new_var_B;
					}
				}
		  
			} // end foreach tmp_num_points outer loop
		
			if (weightsPtr == nullptr) {
				new_mean.alpha *= data_weight;
				new_mean.L *= data_weight;
				new_mean.A *= data_weight;
				new_mean.B *= data_weight;
		  
				new_weight = new_weight_count * data_weight;
		  
				if (!apply_lkm) {
					new_var.alpha *= data_weight;
					new_var.L *= data_weight;
					new_var.A *= data_weight;
					new_var.B *= data_weight;
				}
			}
		
			// Calculate the weight of the old cluster
			double old_weight = total_weight - new_weight; /* weight of C1 */
		
			// Calculate the mean of the new cluster
			new_mean.alpha /= new_weight;
			new_mean.L /= new_weight;
			new_mean.A /= new_weight;
			new_mean.B /= new_weight;    
		
			/* Calculate the mean of the old cluster using the 'combined mean' formula */
			auto& old_mean = mean[old_index]; /* componentwise mean of C1 */
			old_mean.alpha = (total_weight * total_mean.alpha - new_weight * new_mean.alpha) / old_weight;
			old_mean.L = (total_weight * total_mean.L - new_weight * new_mean.L) / old_weight;
			old_mean.A = (total_weight * total_mean.A - new_weight * new_mean.A) / old_weight;
			old_mean.B = (total_weight * total_mean.B - new_weight * new_mean.B) / old_weight;    
		
			/* LOCAL K-MEANS BEGIN */  
			for (int it = 0; it < max_iters; ++it) {
				// Precalculations
				lhs = 0.5 * (sqr(old_mean.alpha) - sqr(new_mean.alpha) + sqr(old_mean.L) - sqr(new_mean.L) + sqr(old_mean.A) - sqr(new_mean.A) + sqr(old_mean.B) - sqr(new_mean.B));
		  
				double rhs_alpha = old_mean.alpha - new_mean.alpha;
				double rhs_L = old_mean.L - new_mean.L;
				double rhs_A = old_mean.A - new_mean.A;
				double rhs_B = old_mean.B - new_mean.B;
		  
				// Reset the statistics of the new cluster
				new_weight = 0.0;
				new_size = 0;
				new_mean.alpha = new_mean.L = new_mean.A = new_mean.B = 0.0;
				new_var.alpha = new_var.L = new_var.A = new_var.B = 0.0;
		  
				for (int ip = 0; ip < tmp_num_points; ) {
					int maxLoopOffset = 0xFFFF;
					int numLeft = (tmp_num_points - ip);
					if (numLeft < maxLoopOffset)
						maxLoopOffset = numLeft;

					maxLoopOffset += ip;
			
					double new_mean_alpha = 0, new_mean_L = 0, new_mean_A = 0, new_mean_B = 0;
			
					double new_var_alpha = 0, new_var_L = 0, new_var_A = 0, new_var_B = 0;
			
					for (; ip < maxLoopOffset; ++ip) {          
						Color c(tmp_data[ip]);
						CIELABConvertor::Lab lab1;
						getLab(c, lab1);
						
						int pointindex = ip;
						if (point_index.get())
							pointindex = point_index[ip];

						if (weightsPtr != nullptr)
							tmp_weight = weightsPtr[pointindex];
			  
						if (lhs < ((rhs_alpha * c.GetA()) + (rhs_L * lab1.L) + (rhs_A * lab1.A) + (rhs_B * lab1.B))) {            
							if (it == max_iters_m1) {
								// Save the membership of the point
								member[pointindex] = old_index;
							}
						}
						else {           
							if (it != max_iters_m1) {
								// Update only mean				  
								if(weightsPtr == nullptr) {
									new_mean_alpha += c.GetA();
									new_mean_L += lab1.L;
									new_mean_A += lab1.A;
									new_mean_B += lab1.B;
								}
								else {
									new_mean.alpha += tmp_weight * c.GetA();
									new_mean.L += tmp_weight * lab1.L;
									new_mean.A += tmp_weight * lab1.A;
									new_mean.B += tmp_weight * lab1.B;
								}
							}
							else {
								// Update mean and variance              
								if (weightsPtr == nullptr) {
									new_mean_alpha += c.GetA();
									new_mean_L += lab1.L;
									new_mean_A += lab1.A;
									new_mean_B += lab1.B;
									
									new_var_alpha += sqr(c.GetA());
									new_var_L += sqr(lab1.L);
									new_var_A += sqr(lab1.A);
									new_var_B += sqr(lab1.B);
								}
								else {
									new_mean.alpha += tmp_weight * c.GetA();
									new_mean.L += tmp_weight * lab1.L;
									new_mean.A += tmp_weight * lab1.A;
									new_mean.B += tmp_weight * lab1.B;
									
									new_var.alpha += tmp_weight * sqr(c.GetA());
									new_var.L += tmp_weight * sqr(lab1.L);
									new_var.A += tmp_weight * sqr(lab1.A);
									new_var.B += tmp_weight * sqr(lab1.B);
								}
					  
								// Save the membership of the point
								member[pointindex] = new_index;
							}
						
							// Update the weight/size of the new cluster					
							if (weightsPtr != nullptr)
								new_weight += tmp_weight;

							++new_size;
						}
					} // end foreach tmp_num_points inner loop
					
					if (weightsPtr == nullptr) {
						new_mean.alpha += new_mean_alpha;
						new_mean.L += new_mean_L;
						new_mean.A += new_mean_A;
						new_mean.B += new_mean_B;
					  
						new_var.alpha += new_var_alpha;
						new_var.L += new_var_L;
						new_var.A += new_var_A;
						new_var.B += new_var_B;
					}        
				} // end foreach tmp_num_points outer loop
				  
				if (weightsPtr == nullptr) {
					new_mean.alpha *= data_weight;
					new_mean.L *= data_weight;
					new_mean.A *= data_weight;
					new_mean.B *= data_weight;
					
					new_weight = new_size * data_weight;
					
					new_var.alpha *= data_weight;
					new_var.L *= data_weight;
					new_var.A *= data_weight;
					new_var.B *= data_weight;
				}
			  
				// Calculate the mean of the new cluster
				new_mean.alpha /= new_weight;
				new_mean.L /= new_weight;
				new_mean.A /= new_weight;
				new_mean.B /= new_weight;
			  
				// Calculate the weight of the old cluster
				old_weight = total_weight - new_weight;
			  
				// Calculate the mean of the old cluster using the 'combined mean' formula
				old_mean.alpha = (total_weight * total_mean.alpha - new_weight * new_mean.alpha) / old_weight;
				old_mean.L = (total_weight * total_mean.L - new_weight * new_mean.L) / old_weight;
				old_mean.A = (total_weight * total_mean.A - new_weight * new_mean.A) / old_weight;
				old_mean.B = (total_weight * total_mean.B - new_weight * new_mean.B) / old_weight;
			}
		
			/* LOCAL K-MEANS END */
		
			/* Store the updated cluster sizes */
			size[old_index] = tmp_num_points - new_size;
			size[new_index] = new_size;
		
			if (new_index == num_colors - 1) {      
				/* This is the last iteration. So, there is no need to determine the cluster to be split in the next iteration. */
				break;
			}
		
			/* Calculate the variance of the new cluster */
			/* Alternative weighted variance formula: ( sum{w_i * x_i^2} / sum{w_i} ) - bar{x}^2 */
			new_var.alpha = new_var.alpha / new_weight - sqr(new_mean.alpha);
			new_var.L = new_var.L / new_weight - sqr(new_mean.L);
			new_var.A = new_var.A / new_weight - sqr(new_mean.A);
			new_var.B = new_var.B / new_weight - sqr(new_mean.B);
		
			/* Calculate the variance of the old cluster using the 'combined variance' formula */
			auto& old_var = var[old_index];
			old_var.alpha = ((total_weight * total_var.alpha -
				new_weight * (new_var.alpha + sqr(new_mean.alpha - total_mean.alpha))) / old_weight) -
				sqr(old_mean.alpha - total_mean.alpha);
				
			old_var.L = ((total_weight * total_var.L -
				new_weight * (new_var.L + sqr(new_mean.L - total_mean.L))) / old_weight) -
				sqr(old_mean.L - total_mean.L);
		
			old_var.A = ((total_weight * total_var.A -
				new_weight * (new_var.A + sqr(new_mean.A - total_mean.A))) / old_weight) -
				sqr(old_mean.A - total_mean.A);
		
			old_var.B = ((total_weight * total_var.B -
				new_weight * (new_var.B + sqr(new_mean.B - total_mean.B))) / old_weight) -
				sqr(old_mean.B - total_mean.B);
		
			/* Store the updated cluster weights */
			weight[old_index] = old_weight;
			weight[new_index] = new_weight;
		
			/* Store the cluster TSEs */
			tse[old_index] = old_weight * (old_var.alpha + old_var.L + old_var.A + old_var.B);
			tse[new_index] = new_weight * (new_var.alpha + new_var.L + new_var.A + new_var.B);
		
			/* STEP 4: DETERMINE THE NEXT CLUSTER TO BE SPLIT */
		
			/* Split the cluster with the maximum TSE */
			max_val = DBL_MIN;
			for (UINT ic = 0; ic <= new_index; ++ic) {
				if (max_val < tse[ic]) {
					max_val = tse[ic];
					old_index = ic;
				}
			}
		
			tmp_num_points = size[old_index];
		
			// Allocate tmp_data and point_index only after initial division and then reuse buffers
		
			if (tmp_buffer_used == 0) {
				// When the initial input points are first split into 2 clusters, allocate tmp_data
				// as a buffer large enough to hold the largest of the 2 initial clusters. This
				// buffer is significantly smaller than the original input size and it can be
				// reused for all smaller cluster sizes.
				
				int largerSize = size[0];
				if (num_colors > 1 && size[1] > largerSize)
					largerSize = size[1];
		  
				tmp_data = tmp_buffer;
		  
				tmp_buffer_used = largerSize;
		  
				// alloc and init to zero
				point_index = make_unique<int[]>(largerSize);
			}
		
			// Setup the points and their indexes in the next cluster to be split    
			int count = 0;			
		
			// Read 1 to N values from member array one at a time
			for (int ip = 0; ip < num_points; ++ip) {
				if (member[ip] == old_index) {   
					tmp_data[count] = data[ip];
					point_index[count++] = ip;
				}
			}
			
			if (count != tmp_num_points) {
				cerr << "Cluster to be split is expected to be of size " << tmp_num_points << " not " << count << " !" << endl;
				return;
			}
		}
	  
		/* Determine the final cluster centers */
		int shift_amount = 8 - num_bits;
		int num_empty = 0; /* # empty clusters */
		UINT colortableOffset = 0;
		for (UINT ic = 0; ic < num_colors; ++ic) {
			if (size[ic] > 0) {
				CIELABConvertor::Lab lab1;
				lab1.alpha = rint(mean[ic].alpha);
				lab1.L = mean[ic].L, lab1.A = mean[ic].A, lab1.B = mean[ic].B;
				Color c(CIELABConvertor::LAB2RGB(lab1));
				
				auto A = ((byte) (mean[ic].alpha + 0.5)) << shift_amount; /* round */
				auto R = ((byte) (c.GetR() + 0.5)) << shift_amount; /* round */
				auto G = ((byte) (c.GetG() + 0.5)) << shift_amount; /* round */
				auto B = ((byte) (c.GetB() + 0.5)) << shift_amount; /* round */
				pPalette->Entries[colortableOffset++] = Color::MakeARGB(A, R, G, B);
			}
			else {
				/* Empty cluster */
				++num_empty;
			}
		}
	  
		if (num_empty)
			cerr << "# empty clusters: " << num_empty << endl;
	  
		nMaxColors = num_colors - num_empty;
	}
	
	static inline bool validate_num_bits(const byte num_bits)
	{
		return (0 < num_bits && num_bits <= 8);
	}
	
	/* TODO: What if num_bits == 0 */

	// This method will reduce the precision of each component of each pixel by setting
	// the number of bits on the right side of the value to zero. Note that this method
	// works properly when inPixels and outPixels are the same buffer to support in
	// place processing.
	void cut_bits(const ARGB* inPixels, const UINT numPixels, ARGB* outPixels,
		const byte num_bits_alpha, const byte num_bits_red, const byte num_bits_green, const byte num_bits_blue)
	{  
		if (!validate_num_bits(num_bits_alpha) || !validate_num_bits(num_bits_red) ||
			!validate_num_bits(num_bits_green) || !validate_num_bits(num_bits_blue))
			return;
  
		byte shift_alpha = 8 - num_bits_alpha;
		byte shift_red = 8 - num_bits_red;
		byte shift_green = 8 - num_bits_green;
		byte shift_blue = 8 - num_bits_blue;
  
		if (shift_alpha == shift_red && shift_alpha == shift_green && shift_alpha == shift_blue) {
			// Shift and mask pixels as whole words when the shift amount
			// for all 4 channel is the same.    
			const UINT shift = shift_red;    
			for (UINT i = 0; i < numPixels; ++i) {
				Color c(inPixels[i]);
				outPixels[i] = Color::MakeARGB(c.GetA() >> shift, c.GetR() >> shift, c.GetG() >> shift, c.GetB() >> shift);
			}
		}
		else {
			for (UINT i = 0; i < numPixels; ++i) {
				Color c(inPixels[i]);
				outPixels[i] = Color::MakeARGB(c.GetA() >> shift_alpha, c.GetR() >> shift_red, c.GetG() >> shift_green, c.GetB() >> shift_blue);
			}
		}
	}

	void DivQuantizer::quant_varpart_fast(const ARGB* inPixels, const UINT numPixels, ColorPalette* pPalette,
		const UINT numRows, const bool allPixelsUnique,
		const int num_bits, const int dec_factor, const int max_iters)
	{	  
		const UINT numCols = numPixels / numRows;
		UINT nMaxColors = pPalette->Count;

		auto inputPixels = make_unique<ARGB[]>(numPixels);
		auto tmpPixels = make_unique<ARGB[]>(numPixels);
	  
		double weightUniform = 0.0;
		unique_ptr<double[]> weightsPtr;
	  
		if (allPixelsUnique && num_bits == 8 && dec_factor == 1) {
			// No duplicate pixels and no decimation or bit shifting
			weightUniform = get_double_scale(numPixels);
			std::copy(inPixels, inPixels + numPixels, inputPixels.get());
		}
		else if (!allPixelsUnique && num_bits == 8) {
			// No cut bits, but duplicate pixels, dedup now
			weightsPtr = calc_color_table(inPixels, numPixels, tmpPixels.get(), numRows, numCols, dec_factor, nMaxColors);
			std::copy(tmpPixels.get(), tmpPixels.get() + numPixels, inputPixels.get());
		}
		else {
			// cut bits with right shift and dedup to generate significantly smaller sized buffer
			cut_bits(inPixels, numPixels, tmpPixels.get(), num_bits, num_bits, num_bits, num_bits);
			weightsPtr = calc_color_table(tmpPixels.get(), numPixels, tmpPixels.get(), numRows, numCols, dec_factor, nMaxColors);
			std::copy(tmpPixels.get(), tmpPixels.get() + numPixels, inputPixels.get());
		}
	  
		if (nMaxColors <= 256)
			DivQuantCluster<byte>(numPixels, inputPixels.get(), tmpPixels.get(), weightUniform, weightsPtr.get(), num_bits, max_iters, pPalette, nMaxColors);
		else
			DivQuantCluster<UINT>(numPixels, inputPixels.get(), tmpPixels.get(), weightUniform, weightsPtr.get(), num_bits, max_iters, pPalette, nMaxColors);
	}
	
	short nearestColorIndex(const ColorPalette* pPalette, const UINT nMaxColors, const ARGB argb)
	{
		short k = 0;
		Color c(argb);

		double mindist = SHORT_MAX;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = 0; i < nMaxColors; i++) {
			Color c2(pPalette->Entries[i]);
			getLab(c2, lab2);

			double curdist = sqr(c2.GetA() - c.GetA());
			if (curdist > mindist)
				continue;

			if (nMaxColors < 256) {
				double deltaL_prime_div_k_L_S_L = CIELABConvertor::L_prime_div_k_L_S_L(lab1, lab2);
				curdist += sqr(deltaL_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double a1Prime, a2Prime, CPrime1, CPrime2;
				double deltaC_prime_div_k_L_S_L = CIELABConvertor::C_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2);
				curdist += sqr(deltaC_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				double barCPrime, barhPrime;
				double deltaH_prime_div_k_L_S_L = CIELABConvertor::H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, barCPrime, barhPrime);
				curdist += sqr(deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;

				curdist += CIELABConvertor::R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
				if (curdist > mindist)
					continue;
			}
			else {
				curdist += sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.B - lab1.B);
				if (curdist > mindist)
					continue;
			}

			mindist = curdist;
			k = i;
		}
		return k;
	}

	bool DivQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT nMaxColors, bool dither)
	{
		const UINT bitmapWidth = pSource->GetWidth();
		const UINT bitmapHeight = pSource->GetHeight();

		int pixelIndex = 0;
		vector<ARGB> pixels(bitmapWidth * bitmapHeight);
		GrabPixels(pSource, pixels, hasSemiTransparency, m_transparentPixelIndex, m_transparentColor);

		auto qPixels = make_unique<short[]>(pixels.size());
		if (nMaxColors > 256) {
			hasSemiTransparency = false;
			dither_image(pixels.data(), nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, qPixels.get(), bitmapWidth, bitmapHeight);
			return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex);
		}

		auto pPaletteBytes = make_unique<byte[]>(pDest->GetPaletteSize());
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;

		if (nMaxColors > 2)
			quant_varpart_fast(pixels.data(), pixels.size(), pPalette);			
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

		if (dither)
			dither_image(pixels.data(), pPalette, nearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight);
		else
			map_colors_mps(pixels.data(), pixels.size(), qPixels.get(), pPalette);
		
		if (m_transparentPixelIndex >= 0) {
			UINT k = qPixels[m_transparentPixelIndex];
			if(nMaxColors > 2)
				pPalette->Entries[k] = m_transparentColor;
			else if (pPalette->Entries[k] != m_transparentColor)
				swap(pPalette->Entries[0], pPalette->Entries[1]);
		}

		return ProcessImagePixels(pDest, pPalette, qPixels.get());
	}

}