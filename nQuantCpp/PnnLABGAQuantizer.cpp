/* Fast pairwise nearest neighbor based genetic algorithm with CIELAB color space genetic algorithm
* Copyright (c) 2023 Miller Cy Chan */

#include "stdafx.h"
#include "PnnLABGAQuantizer.h"
#include "CIELABConvertor.h"
#include "BlueNoise.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <iomanip>
#include <numeric>
#include <random>
#include <mutex>
#include <sstream>
#include <unordered_map>

namespace PnnLABQuant
{	
	UINT _bitmapWidth, _dp = 1;
	UINT _nMaxColors = 256;
	double minRatio = 0, maxRatio = 1.0;

	static unordered_map<string, vector<double> > _fitnessMap;
	static mutex _mutex;

	PnnLABGAQuantizer::PnnLABGAQuantizer(PnnLABQuantizer& pq, Bitmap* pSource, UINT nMaxColors) {
		// increment value when criteria violation occurs
		_objectives.resize(4);
		_bitmapWidth = pSource->GetWidth();
		srand(_bitmapWidth * pSource->GetHeight());
		
		m_pq = make_unique<PnnLABQuantizer>(pq);
		if(pq.IsGA())
			return;

		clear();
		_nMaxColors = nMaxColors;

		bool hasSemiTransparency = false;
		m_pixels.resize(_bitmapWidth * pSource->GetHeight()); 
		m_pq->grabPixels(pSource, m_pixels, _nMaxColors, hasSemiTransparency);
		minRatio = (hasSemiTransparency || nMaxColors < 64) ? .0111 : .85;
		maxRatio = min(1.0, nMaxColors / ((nMaxColors < 64) ? 500.0 : 50.0));
		_dp = maxRatio < .1 ? 10000 : 100;
	}

	PnnLABGAQuantizer::PnnLABGAQuantizer(PnnLABQuantizer& pq, const vector<ARGB>& pixels, const UINT bitmapWidth, UINT nMaxColors)
	{
		m_pq = make_unique<PnnLABQuantizer>(pq);
		// increment value when criteria violation occurs
		_objectives.resize(4);
		m_pixels = pixels;
		_bitmapWidth = bitmapWidth;
		srand(pixels.size());
		_nMaxColors = nMaxColors;
	}

	string PnnLABGAQuantizer::getRatioKey() const
	{
		ostringstream ss;
		ss << (int)(_ratioX * _dp);
		auto difference = abs(_ratioX - _ratioY);
		if (difference <= 0.0000001)
			return ss.str();

		ss << ";" << (int)(_ratioY * _dp * 100);
		return ss.str();
	}

	auto PnnLABGAQuantizer::findByRatioKey(const string& ratioKey) const
	{
		lock_guard<mutex> lock(_mutex);
		auto got = _fitnessMap.find(ratioKey);
		if (got != _fitnessMap.end())
			return got->second;
		return vector<double>();
	}

	using AmplifyFn = function<double(const bool)>;

	void PnnLABGAQuantizer::calculateError(vector<double>& errors) {
		auto maxError = maxRatio < .1 ? .5 : .0625;
		if (m_pq->hasAlpha())
			maxError = 1;

		auto fitness = 0.0;
		bool tooSmall = false; // any error < exp(1.0) concluded as too small
		for (int i = 0; i < errors.size(); ++i) {
			errors[i] /= maxError * m_pixels.size();
			if(errors[i] < 3)
				tooSmall = true;
		}

		AmplifyFn amplifyFn = [tooSmall](const double val) -> double {
			if (tooSmall)
				return exp(val);
			return log(val);
		};

		for (int i = 0; i < errors.size(); ++i) {
			if (i == 0 && errors[i] > maxError)
				errors[i] *= amplifyFn(errors[i]);
			else if (errors[i] > (2 * maxError))
				errors[i] *= amplifyFn(errors[i]);
			fitness -= errors[i];
		}

		_objectives = errors;
		_fitness = fitness;
	}

	void PnnLABGAQuantizer::calculateFitness() {
		auto ratioKey = getRatioKey();
		auto objectives = findByRatioKey(ratioKey);
		if (!objectives.empty()) {
			_fitness = -1.0 * accumulate(objectives.begin(), objectives.end(), 0);
			_objectives = objectives;
			return;
		}

		_objectives.resize(4);
		m_pq->setRatio(_ratioX, _ratioY);
		
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + _nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = _nMaxColors;
		m_pq->pnnquan(m_pixels, pPalette, _nMaxColors);

		auto errors = _objectives;
		fill(errors.begin(), errors.end(), 0);

		int threshold = maxRatio < .1 ? -64 : -112;
		int pixelIndex = 0;
		for (int x = 0; x < m_pixels.size(); ++x, ++pixelIndex)
		{
			if(BlueNoise::RAW_BLUE_NOISE[pixelIndex & 4095] > threshold)
				continue;

			auto argb = m_pixels[pixelIndex];
			Color c(argb);
			CIELABConvertor::Lab lab1, lab2;
			m_pq->getLab(c, lab1);
			auto qPixelIndex = m_pq->nearestColorIndex(pPalette, argb, pixelIndex);
			Color c2(pPalette->Entries[qPixelIndex]);
			m_pq->getLab(c2, lab2);

			if (m_pq->hasAlpha()) {
				errors[0] += sqr(lab2.L - lab1.L);
				errors[1] += sqr(lab2.A - lab1.A);
				errors[2] += sqr(lab2.B - lab1.B);
				errors[3] += sqr(lab2.alpha - lab1.alpha) / exp(1.5);
			}
			else {
				errors[0] += abs(lab2.L - lab1.L);
				errors[1] += sqrt(sqr(lab2.A - lab1.A) + sqr(lab2.B - lab1.B));
			}
		}
		
		calculateError(errors);
		lock_guard<mutex> lock(_mutex);
		_fitnessMap.insert({ ratioKey, _objectives });
	}
	
	bool PnnLABGAQuantizer::QuantizeImage(Bitmap* pDest, bool dither) {
		m_pq->setRatio(_ratioX, _ratioY);
		auto pPaletteBytes = make_unique<BYTE[]>(sizeof(ColorPalette) + _nMaxColors * sizeof(ARGB));
		auto pPalette = (ColorPalette*)pPaletteBytes.get();
		pPalette->Count = _nMaxColors;

		m_pq->pnnquan(m_pixels, pPalette, _nMaxColors);
		return m_pq->QuantizeImage(m_pixels, _bitmapWidth, pPalette, pDest, _nMaxColors, dither);
	}

	void PnnLABGAQuantizer::clear() {
		lock_guard<mutex> lock(_mutex);
		_fitnessMap.clear();
	}

	double randrange(double min, double max)
	{
		auto f = (double) rand() / RAND_MAX;
		return min + f * (max - min);
	}
	
	void PnnLABGAQuantizer::setRatio(double ratioX, double ratioY)
	{
		auto difference = abs(ratioX - ratioY);
		if (difference <= minRatio)
			ratioY = ratioX;
		_ratioX = min(max(ratioX, minRatio), maxRatio);
		_ratioY = min(max(ratioY, minRatio), maxRatio);
	}

	float PnnLABGAQuantizer::getFitness() {
		return (float) _fitness;
	}

	static double rotateLeft(double u, double v, double delta) {
		auto theta = M_PI * randrange(minRatio, maxRatio) / exp(delta);
		auto result = u * sin(theta) + v * cos(theta);
		if(result <= minRatio || result >= maxRatio)
			result = rotateLeft(u, v, delta + .5);
		return result;
	}
	
	static double rotateRight(double u, double v, double delta) {
		auto theta = M_PI * randrange(minRatio, maxRatio) / exp(delta);
		auto result = u * cos(theta) - v * sin(theta);
		if(result <= minRatio || result >= maxRatio)
			result = rotateRight(u, v, delta + .5);
		return result;
	}

	shared_ptr<PnnLABGAQuantizer> PnnLABGAQuantizer::crossover(const PnnLABGAQuantizer& mother, int numberOfCrossoverPoints, float crossoverProbability)
	{
		auto child = makeNewFromPrototype();
		if ((rand() % 100) <= crossoverProbability)
			return child;

		auto ratioX = rotateRight(_ratioX, mother._ratioY, 0.0);
		auto ratioY = rotateLeft(_ratioY, mother._ratioX, 0.0);
		child->setRatio(ratioX, ratioY);
		child->calculateFitness();
		return child;
	}

	static double boxMuller(double value) {
		auto r1 = randrange(minRatio, maxRatio);
		return sqrt(-2 * log(value)) * cos(2 * M_PI * r1);
	}

	bool PnnLABGAQuantizer::dominates(const PnnLABGAQuantizer* right) {
		bool better = false;
		for (int f = 0; f < getObjectives().size(); ++f) {
			if (getObjectives()[f] > right->getObjectives()[f])
				return false;

			if (getObjectives()[f] < right->getObjectives()[f])
				better = true;
		}
		return better;
	}

	void PnnLABGAQuantizer::mutation(int mutationSize, float mutationProbability) {
		// check probability of mutation operation
		if ((rand() % 100) > mutationProbability)
			return;

		auto ratioX = _ratioX;
		auto ratioY = _ratioY;
		if(randrange(.0, 1.0) > .5)
			ratioX = boxMuller(ratioX);
		else
			ratioY = boxMuller(ratioY);

		setRatio(ratioX, ratioY);
		calculateFitness();
	}

	vector<double> PnnLABGAQuantizer::getObjectives() const
	{
		return _objectives;
	}

	vector<double>& PnnLABGAQuantizer::getConvertedObjectives()
	{
		return _convertedObjectives;
	}

	void PnnLABGAQuantizer::resizeConvertedObjectives(int numObj) {
		_convertedObjectives.resize(numObj);
	}

	shared_ptr<PnnLABGAQuantizer> PnnLABGAQuantizer::makeNewFromPrototype() {
		auto child = make_shared<PnnLABGAQuantizer>(*m_pq, m_pixels, _bitmapWidth, _nMaxColors);
		auto minRatio2 = 2.0 * minRatio;
		if(minRatio2 > 1)
			minRatio2 = 0;
		auto ratioX = randrange(minRatio, maxRatio);
		auto ratioY = ratioX < minRatio2 ? randrange(minRatio, maxRatio) : ratioX;
		child->setRatio(ratioX, ratioY);
		child->calculateFitness();
		return child;
	}

	UINT PnnLABGAQuantizer::getMaxColors() const {
		return _nMaxColors;
	}

	
	string PnnLABGAQuantizer::getResult() const
	{
		ostringstream ss;
		auto difference = abs(_ratioX - _ratioY);
		if (difference <= 0.0000001)
			ss << std::setprecision(6) << _ratioX;
		else
			ss << std::setprecision(6) << _ratioX << ", " << _ratioY;
		return ss.str();
	}

}
