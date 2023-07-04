/*
 * Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm Using Reference Point-Based Nondominated Sorting Approach,
 * Part I: Solving Problems With Box Constraints[J]. IEEE Transactions on Evolutionary Computation, 2014, 18(4):577-601.
 * Copyright (c) 2023 Miller Cy Chan
 */

#include "stdafx.h"
#include "NsgaIII.h"
#include "PnnLABGAQuantizer.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <omp.h>
#include <time.h>

using namespace std;

namespace nQuantGA
{
	struct ReferencePoint {
	private:
		int memberSize;
		vector<double> position;
		unordered_map<int, double> potentialMembers;

	public:
		ReferencePoint(int M) {
			memberSize = 0;
			position.resize(M);
			potentialMembers.clear();
		}
		
		inline double& operator[](int index)
		{
			return position[index];
		}

		inline const double& operator[](int index) const
		{
			return position[index];
		}

		void addMember()
		{
			++memberSize;
		}

		void addPotentialMember(int memberInd, double distance)
		{
			auto got = potentialMembers.find(memberInd);
			if (got == potentialMembers.end() || distance < got->second)
				potentialMembers[memberInd] = distance;
		}

		int findClosestMember() const
		{
			auto it = min_element(potentialMembers.begin(), potentialMembers.end(),
				[](const auto& l, const auto& r) { return l.second < r.second; });
			return it->first;
		}

		bool hasPotentialMember() const
		{
			return !potentialMembers.empty();
		}

		int randomMember() const
		{
			if (potentialMembers.empty())
				return -1;

			vector<int> members;
			members.reserve(potentialMembers.size());
			for (const auto& [key, _] : potentialMembers) {
				members.emplace_back(key);
			}
			return members[rand() % potentialMembers.size()];
		}

		void removePotentialMember(int memberInd)
		{
			potentialMembers.erase(memberInd);
		}

		int getMemberSize() const
		{
			return memberSize;
		}
		
		int size() const
		{
			return position.size();
		}

		static void generateRecursive(vector<ReferencePoint>& rps, ReferencePoint& pt, int numObjs, int left, int total, int element) {
			if (element == numObjs - 1) {
				pt[element] = left * 1.0 / total;
				rps.emplace_back(pt);
			}
			else {
				for (int i = 0; i <= left; ++i) {
					pt[element] = i * 1.0 / total;
					generateRecursive(rps, pt, numObjs, left - i, total, element + 1);
				}
			}
		}

		static void generateReferencePoints(vector<ReferencePoint>& rps, int M, const vector<int>& p) {
			ReferencePoint pt(M);
			generateRecursive(rps, pt, M, p[0], p[0], 0);

			if (p.size() > 1) { // two layers of reference points (Check Fig. 4 in NSGA-III paper)
				vector<ReferencePoint> insideRps;
				generateRecursive(insideRps, pt, M, p[1], p[1], 0);

				double center = 1.0 / M;

				for (auto& insideRp : insideRps) {
					for (int j = 0; j < insideRp.getMemberSize(); ++j)
						insideRp[j] = center + insideRp[j] / 2; // (k=num_divisions/M, k, k, ..., k) is the center point

					rps.emplace_back(insideRp);
				}
			}
		}

		static double perpendicularDistance(const ReferencePoint& rp, const vector<double>& point)
		{
			double numerator = 0, denominator = 0;
			for (int i = 0; i < rp.size(); ++i) {
				numerator += rp[i] * point[i];
				denominator += pow(rp[i], 2);
			}

			if (denominator <= 0)
				return (numeric_limits<double>::max)();

			double k = numerator / denominator;
			double d = 0;
			for (int i = 0; i < rp.size(); ++i)
				d += pow(k * rp[i] - point[i], 2);

			return sqrt(d);
		}

	};


	template <class T>
	void associate(vector<ReferencePoint>& rps, const vector<shared_ptr<T> >& pop, const vector<vector<int> >& fronts) {
		for (int t = 0; t < fronts.size(); ++t) {
			for (auto memberInd : fronts[t]) {
				int minRp = rps.size() - 1;
				auto minDist = (numeric_limits<double>::max)();
				for (int r = 0; r < rps.size(); ++r) {
					auto d = ReferencePoint::perpendicularDistance(rps[r], pop[memberInd]->getConvertedObjectives());
					if (d < minDist) {
						minDist = d;
						minRp = r;
					}
				}

				if (t + 1 != fronts.size()) // associating members in St/Fl (only counting)
					rps[minRp].addMember();
				else
					rps[minRp].addPotentialMember(memberInd, minDist);

			}// for - members in front
		}// for - fronts
	}
		

	static unique_ptr<double[]> guassianElimination(vector<vector<double> >& A, const double* b)
	{
		const int N = A.size();
		for (int i = 0; i < N; ++i)
			A[i].emplace_back(b[i]);

		for (int base = 0; base < N - 1; ++base) {
			for (int target = base + 1; target < N; ++target) {
				double ratio = A[target][base] / A[base][base];
				for (int term = 0; term < A[base].size(); ++term)
					A[target][term] -= A[base][term] * ratio;
			}
		}

		auto x = make_unique<double[]>(N);
		for (int i = N - 1; i >= 0; --i) {
			for (int known = i + 1; known < N; ++known)
				A[i][N] -= A[i][known] * x[known];

			x[i] = A[i][N] / A[i][i];
		}
		return x;
	}

	// ----------------------------------------------------------------------
	// ASF: Achivement Scalarization Function
	// ----------------------------------------------------------------------
	static double ASF(const vector<double>& objs, const double* weight)
	{
		auto max_ratio = -(numeric_limits<double>::max)();
		for (int f = 0; f < objs.size(); ++f) {
			auto w = max(weight[f], 1e-6);
			max_ratio = max(max_ratio, objs[f] / w);
		}
		return max_ratio;
	}

	template <class T>
	vector<int> findExtremePoints(const vector<shared_ptr<T> >& pop, const vector<vector<int> >& fronts) {
		const int numObj = pop[0]->getObjectives().size();

		vector<int> exp;
		for (int f = 0; f < numObj; ++f) {
			vector<double> w(numObj, 1e-6);
			w[f] = 1.0;

			auto minASF = (numeric_limits<double>::max)();
			int minIndv = fronts[0].size();

			for (int frontIndv : fronts[0]) { // only consider the individuals in the first front
				auto asf = ASF(pop[frontIndv]->getConvertedObjectives(), w.data());

				if (asf < minASF) {
					minASF = asf;
					minIndv = frontIndv;
				}
			}

			exp.emplace_back(minIndv);
		}

		return exp;
	}

	template <class T>
	vector<double> findMaxObjectives(const vector<shared_ptr<T> >& pop)
	{
		const int numObj = pop[0]->getObjectives().size();
		vector<double> maxPoint(numObj, -(numeric_limits<double>::max)());
		for (int i = 0; i < pop.size(); ++i) {
			for (int f = 0; f < maxPoint.size(); ++f)
				maxPoint[f] = max(maxPoint[f], pop[i]->getObjectives()[f]);
		}

		return maxPoint;
	}

	int findNicheReferencePoint(const vector<ReferencePoint>& rps)
	{
		// find the minimal cluster size
		int minSize = (numeric_limits<int>::max)();
		for (auto& rp : rps)
			minSize = min(minSize, rp.size());

		// find the reference points with the minimal cluster size Jmin
		vector<int> minRps;
		for (int r = 0; r < rps.size(); ++r) {
			if (rps[r].size() == minSize)
				minRps.emplace_back(r);
		}

		// return a random reference point (j-bar)
		return minRps[rand() % minRps.size()];
	}

	template <class T>
	vector<double> constructHyperplane(const vector<shared_ptr<T> >& pop, const vector<int>& extremePoints)
	{
		const int numObj = pop[0]->getObjectives().size();
		// Check whether there are duplicate extreme points.
		// This might happen but the original paper does not mention how to deal with it.
		bool duplicate = false;
		for (int i = 0; !duplicate && i < extremePoints.size(); ++i) {
			for (int j = i + 1; !duplicate && j < extremePoints.size(); ++j)
				duplicate = (extremePoints[i] == extremePoints[j]);
		}

		vector<double> intercepts;

		bool negativeIntercept = false;
		if (!duplicate) {
			// Find the equation of the hyperplane
			vector<double> b(numObj, 1.0);
			vector<vector<double> > A(extremePoints.size());
			for (int p = 0; p < extremePoints.size(); ++p)
				A[p] = pop[ extremePoints[p] ]->getConvertedObjectives();
				
			auto x = guassianElimination(A, b.data());
			// Find intercepts
			for (int f = 0; f < numObj; ++f) {
				intercepts.emplace_back(1.0 / x[f]);

				if(x[f] < 0) {
					negativeIntercept = true;
					break;
				}
			}
		}

		if (duplicate || negativeIntercept) // follow the method in Yuan et al. (GECCO 2015)
			intercepts = findMaxObjectives(pop);
		return intercepts;
	}

	template <class T>
	void normalizeObjectives(vector<shared_ptr<T> >& pop, const vector<vector<int> >& fronts, const vector<double>& intercepts, const vector<double>& idealPoint)
	{
		for (auto& front : fronts) {
			for (int ind : front) {
				auto& convObjs = pop[ind]->getConvertedObjectives();
				for (int f = 0; f < convObjs.size(); ++f) {
					if (abs(intercepts[f] - idealPoint[f]) > 10e-10) // avoid the divide-by-zero error
						convObjs[f] /= intercepts[f] - idealPoint[f];
					else
						convObjs[f] /= 10e-10;
				}
			}
		}
	}

	template <class T>
	vector<vector<int> > nondominatedSort(vector<shared_ptr<T> >& pop) {
		vector<vector<int> > fronts;
		int numAssignedIndividuals = 0;
		int rank = 1;
		auto indvRanks = make_unique<int[]>(pop.size());

		while (numAssignedIndividuals < pop.size()) {
			vector<int> curFront;

			for (int i = 0; i < pop.size(); ++i) {
				if (indvRanks[i] > 0)
					continue; // already assigned a rank

				bool beDominated = false;
				for (int j = 0; j < curFront.size(); ++j) {
					if (pop[curFront[j]]->dominates(pop[i].get()) ) { // i is dominated
						beDominated = true;
						break;
					}
					else if (pop[i]->dominates(pop[ curFront[j] ].get()) ) // i dominates a member in the current front
						curFront.erase(curFront.begin() + j--);
				}

				if (!beDominated)
					curFront.emplace_back(i);
			}

			for (int front : curFront)
				indvRanks[front] = rank;

			fronts.emplace_back(curFront);
			numAssignedIndividuals += curFront.size();
				
			++rank;
		}

		return fronts;
	}

	int selectClusterMember(const ReferencePoint& rp) {
		if (rp.hasPotentialMember()) {
			if (rp.size() == 0) // currently has no member
				return rp.findClosestMember();

			return rp.randomMember();
		}
		return -1;
	}

	template <class T>
	vector<double> translateObjectives(vector<shared_ptr<T> >& pop, const vector<vector<int> >& fronts)
	{
		vector<double> idealPoint;
		const int numObj = pop[0]->getObjectives().size();
		for (int f = 0; f < numObj; ++f) {
			auto minf = (numeric_limits<double>::max)();
			for (int frontIndv : fronts[0]) // min values must appear in the first front
				minf = min(minf, pop[frontIndv]->getObjectives()[f]);
				
			idealPoint.emplace_back(minf);

			for (auto& front : fronts) {
				for (int ind : front) {
					auto chromosome = pop[ind];
					chromosome->resizeConvertedObjectives(numObj);
					chromosome->getConvertedObjectives()[f] = chromosome->getObjectives()[f] - minf;
				}
			}
		}
			
		return idealPoint;
	}

	template <class T>
	vector<shared_ptr<T> > selection(vector<shared_ptr<T> >& cur, vector<ReferencePoint>& rps, const int populationSize) {
		vector<shared_ptr<T> > next;

		// ---------- Step 4 in Algorithm 1: non-dominated sorting ----------
		auto fronts = nondominatedSort(cur);

		// ---------- Steps 5-7 in Algorithm 1 ----------
		int last = 0, next_size = 0;
		while (next_size < populationSize) {
			next_size += fronts[last++].size();
		}

		fronts.resize(last); // remove useless individuals

		for (int t = 0; t < fronts.size() - 1; ++t) {
			for (int frontIndv : fronts[t])
				next.emplace_back(cur[frontIndv]);
		}

		// ---------- Steps 9-10 in Algorithm 1 ----------
		if (next.size() == populationSize)
			return next;

		// ---------- Step 14 / Algorithm 2 ----------
		auto idealPoint = translateObjectives(cur, fronts);

		auto extremePoints = findExtremePoints(cur, fronts);

		auto intercepts = constructHyperplane(cur, extremePoints);

		normalizeObjectives(cur, fronts, intercepts, idealPoint);

		// ---------- Step 15 / Algorithm 3, Step 16 ----------
		associate(rps, cur, fronts);

		// ---------- Step 17 / Algorithm 4 ----------
		while (next.size() < populationSize) {
			int minRp = findNicheReferencePoint(rps);

			int chosen = selectClusterMember(rps[minRp]);
			if (chosen < 0) // no potential member in Fl, disregard this reference point
				rps.erase(rps.begin() + minRp);
			else {
				rps[minRp].addMember();
				rps[minRp].removePotentialMember(chosen);
				next.emplace_back(cur[chosen]);
			}
		}

		return next;
	}

	// Initializes NsgaIII
	template <class T>
	NsgaIII<T>::NsgaIII(T& prototype, int numberOfChromosomes)
	{
		_prototype = prototype.makeNewFromPrototype();

		// there should be at least 2 chromosomes in population
		if (numberOfChromosomes < 2)
			numberOfChromosomes = 2;
		_populationSize = numberOfChromosomes;
	}

	template <class T>
	NsgaIII<T>::NsgaIII(T& prototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability) : NsgaIII<T>(prototype, 9)
	{
		_criteriaLength = prototype.getObjectives().size();
		_mutationSize = mutationSize;
		_numberOfCrossoverPoints = numberOfCrossoverPoints;
		_crossoverProbability = crossoverProbability;
		_mutationProbability = mutationProbability;

		_objDivision.clear();
		if(_criteriaLength < 8)
			_objDivision.emplace_back(6);
		else {
			_objDivision.emplace_back(3);
			_objDivision.emplace_back(2);
		}
	}


	template <class T>
	vector<shared_ptr<T> > NsgaIII<T>::crossing(vector<shared_ptr<T> >& population)
	{
		vector<shared_ptr<T> > offspring;
		offspring.reserve(_populationSize);
		#pragma omp parallel for
		for (int i = 0; i < _populationSize; ++i) {
			if (i % 2 == 0) {
				int father = rand() % _populationSize, mother = rand() % _populationSize;
				offspring.emplace_back(population[father]->crossover(*(population[mother]), _numberOfCrossoverPoints, _crossoverProbability));
				if((i + 1) < _populationSize)
					offspring.emplace_back(population[mother]->crossover(*(population[father]), _numberOfCrossoverPoints, _crossoverProbability));
			}
		}
		return offspring;
	}

	template <class T>
	vector<shared_ptr<T> > NsgaIII<T>::initialize()
	{
		vector<shared_ptr<T> > result(_populationSize);
		result[0] = _prototype->makeNewFromPrototype();
		// initialize new population with chromosomes randomly built using prototype
		#pragma omp parallel for
		for (int i = 1; i < _populationSize; ++i) {
			result[i] = _prototype->makeNewFromPrototype();
		}
		return result;
	}

	template <class T>
	void NsgaIII<T>::reform()
	{
		srand(time(NULL));
		if(_crossoverProbability < 95)
			_crossoverProbability += 1.0f;
		else if(_mutationProbability < 30)
			_mutationProbability += 1.0f;
	}

	template <class T>
	vector<shared_ptr<T> > NsgaIII<T>::replacement(vector<shared_ptr<T> >& population)
	{
		vector<ReferencePoint> rps;
		ReferencePoint::generateReferencePoints(rps, _criteriaLength, _objDivision);
		return selection(population, rps, _populationSize);
	}

	// Starts and executes algorithm
	template <class T>
	void NsgaIII<T>::run(int maxRepeat, double minFitness)
	{
		if (_prototype.get() == nullptr)
			return;

		vector<shared_ptr<T> > pop[2];
		pop[0] = initialize();

		// Current generation
		int currentGeneration = 0;
		int bestNotEnhance = 0;
		double lastBestFit = 0.0;

		int cur = 0, next = 1;
		for (; ;)
		{
			auto best = getResult();
			if(currentGeneration > 0) {
				ostringstream status;
				status << "\rFitness: " << best->getFitness() << "\t Generation: " << currentGeneration;
				wcout << status.str().c_str();

				// algorithm has reached criteria?
				if (best->getFitness() > minFitness)
					break;

				double difference = abs(best->getFitness() - lastBestFit);
				if (difference <= 0.0000001)
					++bestNotEnhance;
				else {
					lastBestFit = best->getFitness();
					bestNotEnhance = 0;
				}

				if (bestNotEnhance > (maxRepeat / 50))
					reform();
			}

			/******************* crossover *****************/
			auto offspring = crossing(pop[cur]);
				
			/******************* mutation *****************/
			#pragma omp parallel for
			for (int i = 0; i < offspring.size(); ++i) {
				offspring[i]->mutation(_mutationSize, _mutationProbability);
			}

			pop[cur].insert(pop[cur].end(), offspring.begin(), offspring.end());

			/******************* replacement *****************/	
			pop[next] = replacement(pop[cur]);
			_best = pop[next][0]->dominates(pop[cur][0].get()) ? pop[next][0] : pop[cur][0];

			swap(cur, next);
			++currentGeneration;
		}
	}

	// explicit instantiations
	template class NsgaIII<PnnLABQuant::PnnLABGAQuantizer>;
}
