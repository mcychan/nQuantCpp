/*
 * Wu, M.; Yang, D.; Zhou, B.; Yang, Z.; Liu, T.; Li, L.; Wang, Z.; Hu,
 * K. Adaptive Population NSGA-III with Dual Control Strategy for Flexible Job
 * Shop Scheduling Problem with the Consideration of Energy Consumption and Weight. Machines 2021, 9, 344.
 * https://doi.org/10.3390/machines9120344
 * Copyright (c) 2023 Miller Cy Chan
 */

#include "stdafx.h"
#include "APNsgaIII.h"
#include "PnnLABGAQuantizer.h"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace nQuantGA
{
	int _currentGeneration = 0, _max_iterations = 5;
	int _maxRepeat = min(15, _max_iterations / 2);

	// Initializes Adaptive Population NSGA-III with Dual Control Strategy
	template <class T>
	APNsgaIII<T>::APNsgaIII(T& prototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability)
		: NsgaIII<T>(prototype, numberOfCrossoverPoints, mutationSize, crossoverProbability, mutationProbability)
	{
	}

	template <class T>
	APNsgaIII<T>::APNsgaIII(T& prototype) : APNsgaIII<T>(prototype, 2, 2, 80.0f, 3.0f)
	{
	}
	
	template <class T>
	double APNsgaIII<T>::ex(T& chromosome)
	{
		double numerator = 0.0, denominator = 0.0;
		for (int f = 0; f < chromosome.getObjectives().size(); ++f) {
			numerator += chromosome.getObjectives()[f] - this->_best->getObjectives()[f];
			denominator += _worst->getObjectives()[f] - this->_best->getObjectives()[f];
		}
		return (numerator + 1) / (denominator + 1);
	}

	template <class T>
	void APNsgaIII<T>::popDec(vector<shared_ptr<T> >& population)
	{
		int N = population.size();
		if(N <= this->_populationSize)
			return;

		int rank = (int) (.3 * this->_populationSize);

		for(int i = 0; i < N; ++i) {
			auto exValue = ex(*population[i]);

			if(exValue > .5 && i > rank) {
				population.erase(population.begin() + i);
				if(--N <= this->_populationSize)
					break;
			}
		}
	}

	template <class T>
	void APNsgaIII<T>::dualCtrlStrategy(vector<shared_ptr<T> >& population, int bestNotEnhance, int nMax)
	{
		int N = population.size();
		int nTmp = N;
		for(int i = 0; i < nTmp; ++i) {
			auto& chromosome = population[i];
			auto pTumor = chromosome->makeNewFromPrototype();
			pTumor->mutation(this->_mutationSize, this->_mutationProbability);

			_worst = population[population.size() - 1];
			if(pTumor->dominates(chromosome.get() )) {
				chromosome = pTumor;
				if(pTumor->dominates(this->_best.get()))
					this->_best = pTumor;
			}
			else {
				if(bestNotEnhance >= _maxRepeat && N < nMax) {
					++N;
					if(_worst->dominates(pTumor.get())) {
						population.emplace_back(pTumor);
						_worst = pTumor;
					}
					else
						population.insert(population.end() - 1, pTumor);
				}
			}
		}
		popDec(population);
	}

	template <class T>
	vector<shared_ptr<T> > APNsgaIII<T>::replacement(vector<shared_ptr<T> >& population)
	{
		auto result = NsgaIII<T>::replacement(population);
		sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs)
		{
			return lhs->getFitness() > rhs->getFitness();
		});
		return result;
	}
		
	// Starts and executes algorithm
	template <class T>
	void APNsgaIII<T>::run(int maxRepeat, double minFitness)
	{
		if (this->_prototype.get() == nullptr)
			return;

		vector<shared_ptr<T> > pop[2];
		pop[0] = this->initialize();
		int nMax = (int) (1.5 * this->_populationSize);

		int bestNotEnhance = 0;
		double lastBestFit = 0.0;

		int cur = 0, next = 1;
		while(_currentGeneration < _max_iterations)
		{
			auto best = this->getResult();
			if(_currentGeneration > 0) {
				double difference = abs(best->getFitness() - lastBestFit);
				if (difference <= 1e-6)
					++bestNotEnhance;
				else {
					lastBestFit = best->getFitness();
					bestNotEnhance = 0;
				}

				ostringstream status;
				if(bestNotEnhance >= _maxRepeat)
					status << "\rFitness: " << showpoint << best->getFitness() << "\t Generation: " << _currentGeneration << " ...";
				else
					status << "\rFitness: " << showpoint << best->getFitness() << "\t Generation: " << _currentGeneration;
				wcout << status.str().c_str();
					
				if (best->getFitness() > minFitness) 
					break;

				if (bestNotEnhance > (maxRepeat / 50))
					this->reform();
			}

			/******************* crossover *****************/
			auto offspring = this->crossing(pop[cur]);
				
			/******************* mutation *****************/
			#pragma omp parallel for
			for (int i = 0; i < offspring.size(); ++i) {
				offspring[i]->mutation(this->_mutationSize, this->_mutationProbability);
			}

			pop[cur].insert(pop[cur].end(), offspring.begin(), offspring.end());
				
			/******************* replacement *****************/		
			pop[next] = replacement(pop[cur]);
			this->_best = pop[next][0]->dominates(pop[cur][0].get()) ? pop[next][0] : pop[cur][0];

			dualCtrlStrategy(pop[next], bestNotEnhance, nMax);

			swap(cur, next);
			++_currentGeneration;
		}

	}

	// explicit instantiations
	template class APNsgaIII<PnnLABQuant::PnnLABGAQuantizer>;
}