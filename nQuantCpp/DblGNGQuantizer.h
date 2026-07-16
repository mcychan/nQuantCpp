#pragma once
#include <memory>
#include <vector>
#include <unordered_map>
#include "CIELABConvertor.h"
using namespace std;

namespace GrowingNeuralGas
{
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

	class DblGNGQuantizer
	{
		private:
			struct GNGNode;

			struct SharedPtrHash {
				template <typename T>
				size_t operator()(const shared_ptr<T>& ptr) const {
					return hash<T*>()(ptr.get());
				}
			};

			struct GNGNode {
				vector<double> weight;
				double error = 0.0;

				unordered_map<shared_ptr<GNGNode>, int, SharedPtrHash> neighbors;

				GNGNode(const vector<double>& w) : weight(w), error(0.0) {}

				void addNeighbour(const shared_ptr<GNGNode>& nextNode) {
					neighbors[nextNode] = 0;
				}

				shared_ptr<GNGNode> findNeighborByMaxError() {
					if (neighbors.empty())
						return nullptr;
					auto it = max_element(neighbors.begin(), neighbors.end(),
						[](const auto& a, const auto& b) -> bool {
							return a.first->error < b.first->error;
						}
					);
					return it->first;
				}

				void incrementAge() {
					for (auto& [neighbor, age] : neighbors) {
						age += 1;
					}
				}

				bool noNeighbor() const { return neighbors.empty(); }
				void removeNeighbour(const shared_ptr<GNGNode>& nextNode) { neighbors.erase(nextNode); }

				void removeNeighbourByAge(int maxAge) {
					for (auto it = neighbors.begin(); it != neighbors.end(); ) {
						if (it->second > maxAge) {
							it = neighbors.erase(it);
						}
						else {
							++it;
						}
					}
				}

				double distance(const vector<double>& input) const {
					double d = 0.0;
					for (size_t i = 0; i < weight.size(); ++i) {
						double diff = weight[i] - input[i];
						d += diff * diff;
					}
					return d;
				}
			};

			bool isGA = false, hasSemiTransparency = false;
			int startingPoints = 2;
			double learningRate = .002;
			double mDivn = 0.0;
			vector<shared_ptr<GNGNode>> nodes;
			vector<shared_ptr<GNGNode>> samples;
			vector<shared_ptr<GNGNode>> uniqueSamples;
			vector<shared_ptr<GNGNode>> stdDevSamples;
			unordered_map<uint32_t, int> histogram;

			unordered_map<int, vector<unsigned short> > closestMap;
			unordered_map<int, unsigned short> nearestMap;
			
			void insertNewNodeWeighted(unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>& assignments);
			void updateNodeWeightsAdaptive(unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>& assignments,
				double baseLearningRate, double progress);
			void manageGraphTopology(unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>& assignments, int remainingEpochs);
			void initializeDistributedNode(const vector<shared_ptr<GNGNode>>& samples, int noOfStartingPoints);
			shared_ptr<GNGNode> findBestWinner(const vector<double>& sample, const vector<shared_ptr<GNGNode>>& snapshot);
			void trainBatch(vector<shared_ptr<GNGNode>>& samples, vector<shared_ptr<GNGNode>>& uniqueSamples,
				vector<shared_ptr<GNGNode>>& stdDevSamples, int totalEpochs);
			void Inxbuild(ARGB* pPalette, const UINT& nMaxColors);

			unsigned short closestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos);
			bool quantize_image(const vector<ARGB>& pixels, const ARGB* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const UINT frameIndex, const bool dither);

		public:
			int m_transparentPixelIndex = -1;
			const double TRANS_RATE = 1 - (512 + 101) / 768.0;
			vector<float> saliencies;
			unordered_map<ARGB, CIELABConvertor::Lab> pixelMap;

			DblGNGQuantizer();
			DblGNGQuantizer(const DblGNGQuantizer& quantizer);
			void clear();
			void gngquan(const vector<ARGB>& pixels, ARGB* pPalette, UINT& nMaxColors);
			const bool IsGA() const;
			void getLab(const Color& c, CIELABConvertor::Lab& lab1);
			const bool hasAlpha() const;
			unsigned short nearestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos);
			void setParams(double learningRate, int startingPoints);
			void grabPixels(Bitmap* srcImg, vector<ARGB>& pixels, UINT& nMaxColors, bool& hasSemiTransparency);
			bool QuantizeImageByPal(const vector<ARGB>& pixels, const UINT bitmapWidth, const ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, const UINT frameIndex = 0, bool dither = true);
			bool QuantizeImage(const vector<ARGB>& pixels, const UINT bitmapWidth, ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither = true);
			bool QuantizeImage(Bitmap* pSource, Bitmap *pDest, UINT& nMaxColors, bool dither = true);
	};
}