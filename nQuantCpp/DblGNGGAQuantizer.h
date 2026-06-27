#pragma once
#include "DblGNGQuantizer.h"
#include "APNsgaIII.h"

#include <string>

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

	class DblGNGGAQuantizer
	{
	private:
		//Asserts floating point compatibility at compile time
		static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");		

		int _startingPoints = 2;
		double _fitness = -numeric_limits<double>::infinity();
		double _learningRate = .002;
		vector<double> _convertedObjectives;
		vector<double> _objectives;
		vector<vector<ARGB> > m_pixelsList;
		unique_ptr<DblGNGQuantizer> m_dgq;

		void calculateError(vector<double>& errors);
		void calculateFitness();
		string getParamsKey() const;
		auto findByParamsKey(const string& paramsKey) const;
		void clear();

	public:
		DblGNGGAQuantizer(DblGNGQuantizer& dgq, const vector<shared_ptr<Bitmap> >& pSources, UINT nMaxColors);
		DblGNGGAQuantizer(DblGNGQuantizer& dgq, const vector<vector<ARGB> >& pixelsList, const vector<UINT>& bitmapWidths, UINT nMaxColors);

		float getFitness();
		shared_ptr<DblGNGGAQuantizer> crossover(const DblGNGGAQuantizer& mother, int numberOfCrossoverPoints, float crossoverProbability);
		bool dominates(const DblGNGGAQuantizer* right);
		void mutation(int mutationSize, float mutationProbability);
		vector<double> getObjectives() const;
		vector<double>& getConvertedObjectives();
		void resizeConvertedObjectives(int numObj);
		shared_ptr<DblGNGGAQuantizer> makeNewFromPrototype();

		UINT getMaxColors() const;
		string getResult() const;
		bool hasAlpha() const {
			return m_dgq->hasAlpha();
		}
		void setParams(double learningRate, int startingPoints);
		bool QuantizeImage(vector<shared_ptr<Bitmap> >& pBitmaps, bool dither = true);
	};

}
