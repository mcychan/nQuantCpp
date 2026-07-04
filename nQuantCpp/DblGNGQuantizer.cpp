#include "stdafx.h"
#include "DblGNGQuantizer.h"
#include "bitmapUtilities.h"
#include "CIELABConvertor.h"
#include "BlueNoise.h"
#include "GilbertCurve.h"
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace GrowingNeuralGas
{
	/* Distributed Batch Learning Growing Neural Gas algorithm with CIELAB color space
	Copyright (c) 2026 Miller Cy Chan
	 * Siow, C. Z., Saputra, A. A., Obo, T., & Kubota, N. (2024).
	 * Distributed batch learning of growing neural gas for quick and efficient clustering.
	 * Mathematics, 12(12), 1909. */

	double PR = 0.299, PG = 0.587, PB = 0.114, PA = .3333;
	BYTE alphaThreshold = 0xF;	

	int maxNodes = 256;		// number of colours used	
	int epochs = 20, maxAge = 10; // Maximum age in terms of epochs/batches, not pixels
	mt19937 random;

	bool enforcedDither = true;
	ARGB m_transparentColor = Color::Transparent;

	static const float coeffs[3][3] = {
		{0.299f, 0.587f, 0.114f},
		{-0.14713f, -0.28886f, 0.436f},
		{0.615f, -0.51499f, -0.10001f}
	};	

	DblGNGQuantizer::DblGNGQuantizer() : startingPoints(2)
		, learningRate(0.002)
		, mDivn(0.0)
		, isGA(false)
	{
	}

	DblGNGQuantizer::DblGNGQuantizer(const DblGNGQuantizer& quantizer) : startingPoints(quantizer.startingPoints)
		, learningRate(quantizer.learningRate)
		, mDivn(quantizer.mDivn)
		, hasSemiTransparency(quantizer.hasSemiTransparency)
		, m_transparentPixelIndex(quantizer.m_transparentPixelIndex)
		, saliencies(quantizer.saliencies)
		, samples(quantizer.samples)
		, uniqueSamples(quantizer.uniqueSamples)
		, stdDevSamples(quantizer.stdDevSamples)
		, isGA(true)
	{
		pixelMap.insert(quantizer.pixelMap.begin(), quantizer.pixelMap.end());
	}	

	int calculateStartingPoints(int transparentPixelIndex, double learningRate) {
		const auto K = (transparentPixelIndex >= 0 || maxNodes > 32) ? 6.5 : 10.0;
		auto continuousPoints = maxNodes * exp(-K * learningRate);
		auto noOfStartingPoints = static_cast<int>(round(continuousPoints));

		auto minFloor = transparentPixelIndex >= 0 ? 8 : 2;
		auto maxCeiling = max(2, maxNodes / 4); 

		return max(minFloor, min(maxCeiling, noOfStartingPoints));
	}

	void DblGNGQuantizer::insertNewNodeWeighted(unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>& assignments) {
		auto it_q = max_element(nodes.begin(), nodes.end(),
			[&](const auto& a, const auto& b) -> bool {
				auto errA = 0.0;
				auto errB = 0.0;

				auto itA = assignments.find(a);
				if (itA != assignments.end() && !itA->second.empty()) {
					errA = a->error / log1p(itA->second.size());
				}

				auto itB = assignments.find(b);
				if (itB != assignments.end() && !itB->second.empty()) {
					errB = b->error / log1p(itB->second.size());
				}

				return errA < errB;
			}
		);

		if (it_q == nodes.end())
			return;

		auto q = *it_q;

		if (q == nullptr || q->noNeighbor())
			return;

		auto f = q->findNeighborByMaxError();
		if (f == nullptr)
			return;

		auto weightSize = q->weight.size();
		auto newWeight = vector<double>(weightSize);
		for (auto i = 0u; i < weightSize; ++i) {
			newWeight[i] = (q->weight[i] + f->weight[i]) / 2.0;
		}

		auto r = make_shared<GNGNode>(newWeight);

		q->removeNeighbour(f);
		f->removeNeighbour(q);

		q->addNeighbour(r);
		f->addNeighbour(r);
		r->addNeighbour(q);
		r->addNeighbour(f);

		nodes.push_back(r);

		q->error *= 0.5;
		f->error *= 0.5;
		r->error = q->error;
	}

	void DblGNGQuantizer::updateNodeWeightsAdaptive(
		unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>& assignments,
		double baseLearningRate, double progress) {
		// Flatten map entries into indexable vector for uniform OpenMP loop parsing
		auto entries = vector<pair<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>>>(assignments.begin(), assignments.end());

		#pragma omp parallel for schedule(dynamic)
		for (auto i = 0; i < static_cast<int>(entries.size()); ++i) {
			auto& node = entries[i].first;
			const auto& cluster = entries[i].second;

			if (cluster.empty()) continue;

			auto mean = vector<double>(node->weight.size(), 0.0);
			for (const auto& sample : cluster) {
				for (auto j = 0u; j < mean.size(); ++j) {
					mean[j] += sample->weight[j];
				}
			}
			for (auto j = 0u; j < mean.size(); ++j) {
				mean[j] /= cluster.size();
			}

			if (progress < 0.4) {
				auto currentLR = baseLearningRate * (1.0 - progress);
				for (auto j = 0u; j < node->weight.size(); ++j) {
					node->weight[j] += currentLR * (mean[j] - node->weight[j]);
				}
			}
			else {
				auto snapFactor = (progress - 0.4) / 0.6;
				for (auto j = 0u; j < node->weight.size(); ++j) {
					node->weight[j] = node->weight[j] + snapFactor * (mean[j] - node->weight[j]);
				}
			}
		}
	}

	void DblGNGQuantizer::manageGraphTopology(unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, 
		SharedPtrHash>& assignments, int remainingEpochs) {
		const auto MAX_AGE = 20;

		for (auto& node : nodes) {
			node->incrementAge();
		}

		for (const auto& [firstWinner, samplesInCluster] : assignments) {
			if (samplesInCluster.empty()) continue;

			const auto& anchorSample = samplesInCluster[0]->weight;
			auto secondWinner = shared_ptr<GNGNode>(nullptr);
			auto minDistance2 = DBL_MAX;

			for (const auto& potentialSecond : nodes) {
				if (potentialSecond == firstWinner) continue;

				auto dist = potentialSecond->distance(anchorSample);
				if (dist < minDistance2) {
					minDistance2 = dist;
					secondWinner = potentialSecond;
				}
			}

			if (secondWinner != nullptr) {
				firstWinner->addNeighbour(secondWinner);
				secondWinner->addNeighbour(firstWinner);
			}
		}

		for (auto& node : nodes) {
			node->removeNeighbourByAge(MAX_AGE);
		}

		nodes.erase(remove_if(nodes.begin(), nodes.end(), [&](const shared_ptr<GNGNode>& node) -> bool {
			if (node->noNeighbor()) {
				auto it = assignments.find(node);
				return it == assignments.end() || it->second.empty();
			}
			return false;
			}), nodes.end());

		auto missingNodes = maxNodes - static_cast<int>(nodes.size());
		if (missingNodes > 0 && remainingEpochs > 0) {
			auto targetInsertions = static_cast<int>(ceil((double)missingNodes / remainingEpochs));
			for (auto i = 0; i < targetInsertions; ++i) {
				if (static_cast<int>(nodes.size()) < maxNodes) {
					insertNewNodeWeighted(assignments);
				}
				else {
					break;
				}
			}
		}
	}

	void DblGNGQuantizer::getLab(const Color& c, CIELABConvertor::Lab& lab1)
	{
		auto got = pixelMap.find(c.GetValue());
		if (got == pixelMap.end()) {
			CIELABConvertor::RGB2LAB(c, lab1);
			pixelMap[c.GetValue()] = lab1;
		}
		else
			lab1 = got->second;
	}

	void DblGNGQuantizer::initializeDistributedNode(const vector<shared_ptr<GNGNode>>& samples, int noOfStartingPoints) {
		if (samples.empty()) {
			throw invalid_argument("Sample list cannot be empty.");
		}

		auto actualStartingPoints = min(noOfStartingPoints, static_cast<int>(samples.size()));
		if (actualStartingPoints < 2) actualStartingPoints = 2;

		maxAge = max(maxAge, actualStartingPoints * (epochs / 10));
		random.seed(samples.size());
		nodes.clear();

		auto initialNodes = vector<shared_ptr<GNGNode>>();
		auto chosenIndices = unordered_set<int>();
		auto dist = uniform_int_distribution<int>(0, samples.size() - 1);

		while (static_cast<int>(initialNodes.size()) < actualStartingPoints) {
			auto randomIndex = dist(random);
			if (chosenIndices.insert(randomIndex).second) {
				initialNodes.push_back(make_shared<GNGNode>(samples[randomIndex]->weight));
			}
		}

		for (auto i = 0u; i < initialNodes.size(); ++i) {
			auto currentNode = initialNodes[i];
			auto nextNode = initialNodes[(i + 1) % initialNodes.size()];
			currentNode->neighbors[nextNode] = 0;
			nextNode->neighbors[currentNode] = 0;
		}

		nodes.insert(nodes.end(), initialNodes.begin(), initialNodes.end());
	}

	double calculateBalancedLearningRate(int stdDevSampleSize) {
		const auto EPSILON_BASE = 0.12; 
		if (stdDevSampleSize <= 0)
			return EPSILON_BASE;

		auto ratio = (double) maxNodes / stdDevSampleSize;
		auto adaptiveLR = EPSILON_BASE * sqrt(ratio);

		return max(0.015, min(0.080, adaptiveLR));
	}

	shared_ptr<DblGNGQuantizer::GNGNode> DblGNGQuantizer::findBestWinner(const vector<double>& sample, const vector<shared_ptr<GNGNode>>& snapshot) {
		auto winner = shared_ptr<GNGNode>(nullptr);
		auto minDist = DBL_MAX;
		for (auto i = 0u; i < snapshot.size(); i++) {
			auto d = snapshot[i]->distance(sample);
			if (d < minDist) {
				minDist = d;
				winner = snapshot[i];
			}
		}
		if (winner != nullptr) {
			#pragma omp critical(WinnerErrorAccumulation)
			{
				winner->error += minDist;
			}
		}
		return winner;
	}

	void DblGNGQuantizer::trainBatch(vector<shared_ptr<GNGNode>>& samples, vector<shared_ptr<GNGNode>>& uniqueSamples, 
		vector<shared_ptr<GNGNode>>& stdDevSamples, int totalEpochs)
	{
		if (!isGA) {
			learningRate = calculateBalancedLearningRate(stdDevSamples.size());
			startingPoints = calculateStartingPoints(m_transparentPixelIndex, learningRate);
		}

		initializeDistributedNode(uniqueSamples, startingPoints);
		auto growthEpochs = static_cast<int>(totalEpochs * 0.7);

		// PHASE 2: TOPOLOGY GROWTH
		for (auto epoch = 0; epoch < growthEpochs; ++epoch) {
			for (auto& node : nodes) node->error = 0.0;
			auto currentNodesSnapshot = nodes;

			unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash> assignments;
			for (const auto& sample : stdDevSamples) {
				auto winner = findBestWinner(sample->weight, currentNodesSnapshot);
				if (winner) {
					assignments[winner].push_back(sample);
				}
			}

			auto progress = (double)epoch / growthEpochs;
			updateNodeWeightsAdaptive(assignments, learningRate, progress);
			manageGraphTopology(assignments, growthEpochs - epoch);
		}

		// PHASE 3: FINAL CENTROID TUNING
		auto tuningEpochs = totalEpochs - growthEpochs;
		for (auto epoch = 0; epoch < tuningEpochs; ++epoch) {
			auto currentNodesSnapshot = nodes;
			auto realAssignments = unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>();

			// Concurrent collection parallelized securely via OpenMP thread-private map reductions
			#pragma omp parallel
			{
				auto localAssignments = unordered_map<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>, SharedPtrHash>();

				#pragma omp for nowait
				for (auto i = 0; i < static_cast<int>(samples.size()); ++i) {
					auto winner = findBestWinner(samples[i]->weight, currentNodesSnapshot);
					if (winner) localAssignments[winner].push_back(samples[i]);
				}

				#pragma omp critical(MergeRealAssignments)
				{
					for (auto& [node, cluster] : localAssignments) {
						auto& globalCluster = realAssignments[node];
						globalCluster.insert(globalCluster.end(), cluster.begin(), cluster.end());
					}
				}
			}

			// Equivalent to entrySet().parallelStream().forEach(...)
			auto realEntries = vector<pair<shared_ptr<GNGNode>, vector<shared_ptr<GNGNode>>>>(realAssignments.begin(), realAssignments.end());

			#pragma omp parallel for schedule(dynamic)
			for (auto i = 0; i < static_cast<int>(realEntries.size()); ++i) {
				auto& node = realEntries[i].first;
				const auto& cluster = realEntries[i].second;
				if (cluster.empty()) continue;

				auto trueMean = vector<double>(node->weight.size(), 0.0);
				for (const auto& n : cluster) {
					for (auto j = 0u; j < trueMean.size(); ++j) {
						trueMean[j] += n->weight[j];
					}
				}

				for (auto j = 0u; j < node->weight.size(); ++j) {
					node->weight[j] = trueMean[j] / cluster.size();
				}
			}

			// Topology management and aged-link pruning
			for (auto& node : nodes) {
				node->incrementAge();
			}

			for (auto& node : nodes) {
				node->removeNeighbourByAge(maxAge);
			}

			nodes.erase(remove_if(nodes.begin(), nodes.end(), [&](const auto& node) {
				if (node->noNeighbor()) {
					auto it = realAssignments.find(node);
					return it == realAssignments.end() || it->second.empty();
				}
				return false;
				}), nodes.end());

			auto missingNodes = maxNodes - static_cast<int>(nodes.size());
			if (missingNodes > 0 && epochs - epoch > 0) {
				auto targetInsertions = static_cast<int>(ceil((double)missingNodes / (epochs - epoch)));
				for (auto i = 0; i < targetInsertions; ++i) {
					if (static_cast<int>(nodes.size()) >= maxNodes) break;
					insertNewNodeWeighted(realAssignments);
				}
			}
		}
	}

	void DblGNGQuantizer::Inxbuild(ARGB* pPalette, const UINT& nMaxColors) {
		UINT k = 0;

		for (const auto& n : nodes) {
			const auto& channels = n.get()->weight;
			CIELABConvertor::Lab lab1;
			lab1.L = (float) channels[0], lab1.A = (float) channels[1], lab1.B = (float) channels[2];
			if (channels.size() > 3)
				lab1.alpha = (float) channels[3];
			pPalette[k++] = CIELABConvertor::LAB2RGB(lab1);
			if (k >= nMaxColors)
				break;
		}
	}

	unsigned short DblGNGQuantizer::nearestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		int offset = GetARGBIndex(argb, hasSemiTransparency, hasAlpha());
		auto got = nearestMap.find(offset);
		if (got != nearestMap.end())
			return got->second;

		unsigned short k = 0;
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			c = m_transparentColor;

		if (nMaxColors > 2 && hasAlpha() && c.GetA() > alphaThreshold)
			k = 1;

		double mindist = 1e100;
		CIELABConvertor::Lab lab1, lab2;
		getLab(c, lab1);

		for (UINT i = k; i < nMaxColors; ++i) {
			Color c2(pPalette[i]);
			auto curdist = hasSemiTransparency ? sqr(c2.GetA() - c.GetA()) * TRANS_RATE : 0;
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
			else if (hasSemiTransparency || nMaxColors < 16) {
				curdist += sqr(lab2.L - lab1.L);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.A - lab1.A);
				if (curdist > mindist)
					continue;

				curdist += sqr(lab2.B - lab1.B);
			}
			else if (nMaxColors > 32) {
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
		nearestMap[offset] = k;
		return k;
	}
	
	unsigned short DblGNGQuantizer::closestColorIndex(const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos)
	{
		Color c(argb);
		if (c.GetA() <= alphaThreshold)
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);

		vector<unsigned short> closest(4);
		int offset = GetARGBIndex(argb, hasSemiTransparency, hasAlpha());
		auto got = closestMap.find(offset);
		if (got == closestMap.end()) {
			closest[2] = closest[3] = USHRT_MAX;

			for (UINT k = 0; k < nMaxColors; ++k) {
				Color c2(pPalette[k]);

				auto err = PR * sqr(c2.GetR() - c.GetR());
				if (err >= closest[3])
					continue;

				err += PG * sqr(c2.GetG() - c.GetG());
				if (err >= closest[3])
					continue;

				err += PB * sqr(c2.GetB() - c.GetB());
				if (err >= closest[3])
					continue;

				if (hasSemiTransparency)
					err += PA * sqr(c2.GetA() - c.GetA());

				if (err < closest[2]) {
					closest[1] = closest[0];
					closest[3] = closest[2];
					closest[0] = k;
					closest[2] = err;
				}
				else if (err < closest[3]) {
					closest[1] = k;
					closest[3] = err;
				}
			}

			if (closest[3] == USHRT_MAX)
				closest[1] = closest[0];

			closestMap[offset] = closest;
		}
		else
			closest = got->second;

		int idx = 1;
		if (closest[2] == 0 || (rand() % (int)ceil(closest[3] + closest[2])) <= closest[3])
			idx = 0;

		auto MAX_ERR = nMaxColors;
		Color pixel(pPalette[closest[idx]]);
		if (closest[idx + 2] >= MAX_ERR || closest[idx] == 0 || pixel.GetA() < c.GetA())
			return nearestColorIndex(pPalette, nMaxColors, argb, pos);
		return closest[idx];
	}

	void DblGNGQuantizer::gngquan(const vector<ARGB>& pixels, ARGB* pPalette, UINT& nMaxColors)
	{
		maxNodes = nMaxColors;		// number of colours used		

		if (samples.empty()) {
			// Sequential color extraction loop
			for (auto i = 0u; i < pixels.size(); ++i) {
				auto pixel = GetARGB(pixels[i], hasSemiTransparency, hasAlpha());
				auto isRegistered = pixelMap.find(pixel) != pixelMap.end();

				Color c(pixel);
				CIELABConvertor::Lab lab1;
				getLab(c, lab1);
				vector<double> currentWeight;

				if (hasAlpha()) {
					currentWeight = { lab1.L, lab1.A, lab1.B, (double)lab1.alpha };
				}
				else {
					currentWeight = { lab1.L, lab1.A, lab1.B };
				}

				// Generate a new GNGNode inside a safe shared heap block
				auto sampleNode = make_shared<GNGNode>(currentWeight);
				samples.push_back(sampleNode);

				if (!isRegistered) {
					uniqueSamples.push_back(make_shared<GNGNode>(currentWeight));
				}

				// Idiomatic count tracking incrementation (Auto-initializes to 0 if key is absent)
				histogram[pixel]++;
			}

			if (pixelMap.size() <= nMaxColors) {
				/* Fill palette */
				nMaxColors = pixelMap.size();
				int k = 0;
				for (const auto& [pixel, lab] : pixelMap) {
					Color c(pPalette[k]);
					pPalette[k++] = pixel;

					if (k > 1 && c.GetA() == 0)
						swap(pPalette[k - 1], pPalette[0]);
				}

				return;
			}

			mDivn = min(0.9, nMaxColors * 1.0 / pixelMap.size());
			if (hasSemiTransparency)
				mDivn *= -1;
			
			for (const auto& [pixel, count] : histogram) {
				auto freq = static_cast<int>(sqrt(count));

				for (auto j = 0; j < freq; ++j) {
					Color c(pixel);
					CIELABConvertor::Lab lab1;
					getLab(c, lab1);

					vector<double> currentWeight;
					if (hasAlpha()) {
						currentWeight = { lab1.L, lab1.A, lab1.B, (double)lab1.alpha };
					}
					else {
						currentWeight = { lab1.L, lab1.A, lab1.B };
					}

					stdDevSamples.push_back(make_shared<GNGNode>(currentWeight));
				}
			}
		}

		if (mDivn > .0029 && nMaxColors <= 32)
			enforcedDither = false;

		if ((nMaxColors < 32 && mDivn > .015 && mDivn < .032) || (nMaxColors >= 32 && nMaxColors < 64 && mDivn > .03 && mDivn < .06))
			trainBatch(uniqueSamples, samples, stdDevSamples, epochs);
		else
			trainBatch(samples, uniqueSamples, stdDevSamples, epochs);

		if (nodes.size() > static_cast<size_t>(nMaxColors)) {
			wcerr << L"Truncated no. of clusters from " << nodes.size() << L" to " << nMaxColors << L"\n";
		} 
		else if (nodes.size() < static_cast<size_t>(nMaxColors)) {
			nMaxColors = nodes.size();
			wcerr << L"Reduced no. of clusters to " << nMaxColors << L"\n";
		}

		Inxbuild(pPalette, nMaxColors);
	}	

	void DblGNGQuantizer::grabPixels(Bitmap* srcImg, vector<ARGB>& pixels, UINT& nMaxColors, bool& hasSemiTransparency)
	{
		int semiTransCount = 0;
		GrabPixels(srcImg, pixels, semiTransCount, m_transparentPixelIndex, m_transparentColor, alphaThreshold, nMaxColors);
		this->hasSemiTransparency = hasSemiTransparency = semiTransCount > 0;
	}

	bool DblGNGQuantizer::quantize_image(const vector<ARGB>& pixels, const ARGB* pPalette, const UINT nMaxColors, unsigned short* qPixels, const UINT width, const UINT height, const bool dither)
	{
		if (dither) {
			auto NearestColorIndex = [this, nMaxColors](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
				if (nMaxColors <= 4)
					return nearestColorIndex(pPalette, nMaxColors, argb, pos);
				return closestColorIndex(pPalette, nMaxColors, argb, pos);
			};
			return dither_image(pixels.data(), pPalette, nMaxColors, NearestColorIndex, hasSemiTransparency, m_transparentPixelIndex, qPixels, width, height, saliencies, enforcedDither);
		}

		UINT pixelIndex = 0;
		for (UINT j = 0; j < height; ++j) {
			for (UINT i = 0; i < width; ++i)
				qPixels[pixelIndex++] = nearestColorIndex(pPalette, nMaxColors, pixels[pixelIndex], i + j);
		}

		return true;
	}

	bool DblGNGQuantizer::QuantizeImageByPal(const vector<ARGB>& pixels, const UINT bitmapWidth, const ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapHeight = pixels.size() / bitmapWidth;

		if (dither) {
			saliencies.resize(pixels.size());
			auto saliencyBase = .1f;
			for (int i = 0; i < pixels.size(); ++i) {
				const auto& pixel = pixels[i];
				Color c(pixel);

				CIELABConvertor::Lab lab1;
				getLab(c, lab1);

				saliencies[i] = saliencyBase + (1 - saliencyBase) * lab1.L / 100.0f * lab1.alpha / 255.0f;
			}
		}

		if (enforcedDither)
			enforcedDither = nMaxColors < 32 || nMaxColors > 64;

		bool sortedByYDiff = nMaxColors >= 128 && mDivn >= .02 && (!hasAlpha() || mDivn < .18);
		if (dither && (sortedByYDiff || !enforcedDither)) {
			auto qPixels = make_unique<unsigned short[]>(pixels.size());
			quantize_image(pixels, pPalette, nMaxColors, qPixels.get(), bitmapWidth, bitmapHeight, dither);

			pixelMap.clear();
			clear();

			return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex >= 0);
		}		

		auto GetColorIndex = [&](const Color& c) -> int {
			return GetARGBIndex(c, hasSemiTransparency, hasAlpha());
		};
		auto NearestColorIndex = [this, nMaxColors](const ARGB* pPalette, const UINT nMaxColors, ARGB argb, const UINT pos) -> unsigned short {
			if (nMaxColors <= 4)
				return nearestColorIndex(pPalette, nMaxColors, argb, pos);
			return closestColorIndex(pPalette, nMaxColors, argb, pos);
		};

		if (nMaxColors > 256) {
			auto qPixels = make_unique<ARGB[]>(pixels.size());
			Peano::GilbertCurve::dithering(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, NearestColorIndex, GetColorIndex, qPixels.get(), saliencies.data(), mDivn, dither, enforcedDither);

			pixelMap.clear();
			clear();
			return ProcessImagePixels(pDest, qPixels.get(), hasSemiTransparency, m_transparentPixelIndex);
		}

		auto qPixels = make_unique<unsigned short[]>(pixels.size());
		Peano::GilbertCurve::dithering(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, NearestColorIndex, GetColorIndex, qPixels.get(), saliencies.data(), mDivn, dither, enforcedDither);

		if (!dither && nMaxColors > 32) {
			const auto delta = sqr(nMaxColors) / pixelMap.size();
			mDivn = delta > 0.023 ? 1.0f : (float)(36.921 * delta + 0.906);
			BlueNoise::dither(bitmapWidth, bitmapHeight, pixels.data(), pPalette, nMaxColors, NearestColorIndex, GetColorIndex, qPixels.get(), mDivn);
		}

		pixelMap.clear();
		clear();

		return ProcessImagePixels(pDest, qPixels.get(), m_transparentPixelIndex >= 0);
	}

	bool DblGNGQuantizer::QuantizeImage(const vector<ARGB>& pixels, const UINT bitmapWidth, ARGB* pPalette, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		if (nMaxColors <= 32)
			PR = PG = PB = PA = 1;
		else {
			PR = coeffs[0][0]; PG = coeffs[0][1]; PB = coeffs[0][2];
		}

		if (nMaxColors > 2)
			gngquan(pixels, pPalette, nMaxColors);
		else {
			if (m_transparentPixelIndex >= 0) {
				pPalette[0] = m_transparentColor;
				pPalette[1] = Color::Black;
			}
			else {
				pPalette[0] = Color::Black;
				pPalette[1] = Color::White;
			}
		}

		if (m_transparentPixelIndex >= 0) {
			auto k = nearestColorIndex(pPalette, nMaxColors, pixels[m_transparentPixelIndex], m_transparentPixelIndex);
			if (nMaxColors > 2)
				pPalette[k] = m_transparentColor;
			else if (pPalette[k] != m_transparentColor)
				swap(pPalette[0], pPalette[1]);
		}

		const auto& pPal = pPalette;
		return QuantizeImageByPal(pixels, bitmapWidth, pPal, pDest, nMaxColors, dither);
	}

	bool DblGNGQuantizer::QuantizeImage(Bitmap* pSource, Bitmap* pDest, UINT& nMaxColors, bool dither)
	{
		const auto bitmapWidth = pSource->GetWidth();
		const auto bitmapHeight = pSource->GetHeight();
		const auto area = (size_t)(bitmapWidth * bitmapHeight);

		vector<ARGB> pixels(area);
		int semiTransCount = 0;
		grabPixels(pSource, pixels, nMaxColors, hasSemiTransparency);

		if (nMaxColors > 256) {
			auto pPalettes = make_unique<ARGB[]>(nMaxColors);
			auto pPalette = pPalettes.get();
			return QuantizeImage(pixels, bitmapWidth, pPalette, pDest, nMaxColors, dither);
		}

		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = nMaxColors;
		auto result = QuantizeImage(pixels, bitmapWidth, pPalette->Entries, pDest, nMaxColors, dither);
		pDest->SetPalette(pPalette);
		return result;
	}

	void DblGNGQuantizer::clear() {
		saliencies.clear();
		closestMap.clear();
		nearestMap.clear();
	}
	
	const bool DblGNGQuantizer::IsGA() const {
		return isGA;
	}

	const bool DblGNGQuantizer::hasAlpha() const {
		return m_transparentPixelIndex >= 0;
	}

	void DblGNGQuantizer::setParams(double learningRate, int startingPoints) {
		this->learningRate = learningRate;
		this->startingPoints = startingPoints;
		closestMap.clear();
		nearestMap.clear();
	}
}
