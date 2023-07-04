#pragma once
#include "NsgaIII.h"

namespace nQuantGA
{
	/*
	 * Wu, M.; Yang, D.; Zhou, B.; Yang, Z.; Liu, T.; Li, L.; Wang, Z.; Hu,
	 * K. Adaptive Population NSGA-III with Dual Control Strategy for Flexible Job
	 * Shop Scheduling Problem with the Consideration of Energy Consumption and Weight. Machines 2021, 9, 344.
	 * https://doi.org/10.3390/machines9120344
	 * Copyright (c) 2023 Miller Cy Chan
	 */
 
	template <class T>
	class APNsgaIII : public NsgaIII<T>
	{
	protected:
		// Worst of chromosomes
		shared_ptr<T> _worst;

		void dualCtrlStrategy(vector<shared_ptr<T> >& population, int bestNotEnhance, int nMax);
		double ex(T& chromosome);
		void popDec(vector<shared_ptr<T> >& population);
		vector<shared_ptr<T> > replacement(vector<shared_ptr<T> >& population) override;

	public:
		APNsgaIII<T>(T& pPrototype, int numberOfCrossoverPoints, int mutationSize, float crossoverProbability, float mutationProbability);
		APNsgaIII<T>(T& pPrototype);

		// Starts and executes algorithm
		void run(int maxRepeat, double minFitness) override;
	};

}
