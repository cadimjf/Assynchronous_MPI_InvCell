#ifndef OPERATORS_H_
#define OPERATORS_H_
#include <ctime>

#include <iostream>
#include <cstdio>
#include "RealPopulation.h"

class Operators {
	public:
		static int coin();
		static RealGenome* selectParents(RealPopulation &b, int optmizationType);
		static RealGenome* selectParentsByRank(RealPopulation &b);
		static RealGenome* selectParentsByRoullette(RealPopulation &pop, double sum);
		static RealGenome  crossover(RealGenome *parents, double probCrossover);
		static RealGenome* blendCrossover(RealGenome x, RealGenome y, double probCrossover);
		static RealPopulation restorePopulation(string fileName, int nInd, int nP);
		static void mutation(RealGenome &ind, double probMutation, double LB, double UB, int t, int T);
		static void mutation(RealGenome &ind, double probMutation, double *LB, double *UB, int t, int T);
		static double evaluateAllFrequencies(RealGenome ind, char* filename);
		static double evaluate(RealGenome ind);
		static void quickSort(RealPopulation &pop, int low, int sup, int optType );
		static int randI( int min , int max);
		static double randD( double min, double max);
		static void saveFitnessMap(RealPopulation &pop, double *bestParams, int generation, int n, double *l, double *r );
		static void savePopulation(RealPopulation &pop, int gen);
		enum OptmizationType { MAX, MIN };
};

#endif /*OPERATORS_H_*/

