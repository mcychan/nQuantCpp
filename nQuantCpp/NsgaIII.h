#pragma once

#include <memory>
#include <vector>
using namespace std;

namespace nQuantGA
{
	/*
	 * Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm Using Reference Point-Based Nondominated Sorting Approach,
	 * Part I: Solving Problems With Box Constraints[J]. IEEE Transactions on Evolutionary Computation, 2014, 18(4):577-601.
	 * Copyright (c) 2023 Miller Cy Chan
	 */
	template <class T>
	class NsgaIII
	{
	private:		
		// Initializes NsgaIII
		NsgaIII(T& prototype, int numberOfChromosomes);

	protected:
		// Prototype of chromosomes in population
		shared_ptr<T> _prototype;

		// Number of chromosomes
		int _populationSize;

		// Number of crossover points of parent's class tables
		int _numberOfCrossoverPoints;

		// Number of classes that is moved randomly by single mutation operation
		int _mutationSize;

		// Probability that crossover will occur
		float _crossoverProbability;

		// Probability that mutation will occur
		float _mutationProbability;

		vector<int> _objDivision;

		int _criteriaLength;

		shared_ptr<T> _best;

		virtual vector<shared_ptr<T> > crossing(vector<shared_ptr<T> >& population);
		virtual vector<shared_ptr<T> > initialize();
		virtual void reform();
		virtual vector<shared_ptr<T> > replacement(vector<shared_ptr<T> >& population);

	public:
		NsgaIII<T>(T& pPrototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability);

		// Returns pointer to best chromosomes in population
		T* getResult() const
		{
			return _best.get();
		}

		// Starts and executes algorithm
		virtual void run(int maxRepeat, double minFitness);
	};

}
