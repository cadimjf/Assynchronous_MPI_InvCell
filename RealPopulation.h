#ifndef RealPopulation_H_
#define RealPopulation_H_
#include "RealGenome.h"

#include <vector>

class RealPopulation {
public:
	RealPopulation(int numIndividuals, int size, double l, double r);
	RealPopulation(int numIndividuals, int size, double *l, double *r);
	RealPopulation(int numIndividuals, int size);
	RealPopulation(const RealPopulation &m);
	virtual ~RealPopulation();
	void initializeRandomPopulation(double lowerBound, double upperBound);
	void initializeRandomPopulation(double *lowerBound, double *upperBound);
	void initializeZeroPopulation();
	RealGenome &getIndividual(int i);
	int getNumIndividuals();
	void setIndividual(int i, RealGenome ind);
	void copyPopulation(RealPopulation &b);
	RealGenome getBestIndividual(int optimizationType);
	const RealPopulation& operator = (const RealPopulation & rhs);
	void setIndFitness(int index, double value);
	double getSize();
	void evaluateAll(int n, double leftBoundary, double rightBoundary);
	double getFitnessSum();
    int getWorstIndividualIndex(int optType);
    void printPopulation();
	
private:
	vector<RealGenome> population;
	int numIndividuals;
	int size;
	void allocPopulation(int n);
	void copyPop(vector<RealGenome> dest, vector<RealGenome> src);
	
	
};

#endif /*RealPopulation_H_*/
