#include "RealPopulation.h"
#include "Operators.h"

#include <iostream>

using namespace std;

RealPopulation::RealPopulation(int numIndividuals, int size, double l, double r) {
	
		this->numIndividuals = numIndividuals;
		this->size = size;	
		initializeRandomPopulation(l,r);
}

RealPopulation::RealPopulation(int numIndividuals, int size, double *l, double *r) {
	
		this->numIndividuals = numIndividuals;
		this->size = size;	
		initializeRandomPopulation(l,r);
}

RealPopulation::RealPopulation(int numIndividuals, int size) {
	
		this->numIndividuals = numIndividuals;
		this->size = size;	
		initializeZeroPopulation();
}

RealPopulation::RealPopulation(const RealPopulation &m)
{
	numIndividuals = m.numIndividuals;	
	copyPop(population, m.population);
	
}

RealPopulation::~RealPopulation() {	
	
	//delete population;
}

void RealPopulation::initializeRandomPopulation(double lowerBound, double upperBound) {
	
	RealGenome b;

	for(int i = 0; i < numIndividuals; i++) {
			b.setSize(size);
			b.initialize(lowerBound, upperBound);
			population.push_back(b);			
	}

}

void RealPopulation::initializeRandomPopulation(double *lowerBound, double *upperBound) {
	
	RealGenome b;

	for(int i = 0; i < numIndividuals; i++) {
			b.setSize(size);
			b.initialize(lowerBound, upperBound);
			population.push_back(b);			
	}

}

void RealPopulation::initializeZeroPopulation() {
	
	RealGenome b;

	for(int i = 0; i < numIndividuals; i++) {
			b.setSize(size);
			b.initializeZeros();
			population.push_back(b);			
	}

}

RealGenome &RealPopulation::getIndividual(int i) {

	return population.at(i);

}

int RealPopulation::getNumIndividuals() {
	
	return numIndividuals;

}

void RealPopulation::setIndividual(int i, RealGenome ind) {

	population.at(i) = ind;

}

void RealPopulation::copyPop(vector<RealGenome> dest, vector<RealGenome> src) {
	
	dest.clear();

	for(unsigned int i = 0; i < src.size(); i++) {						
			dest.push_back(src.at(i));
	}
}

RealGenome RealPopulation::getBestIndividual(int optimizationType) {


	double best, aux;
	int iBest = 0;
	
	best = population[0].getFitness();

	if(optimizationType == Operators::MAX) {
		for(int i = 1; i < this->getNumIndividuals(); i++) {
			aux = population.at(i).getFitness();
			if(aux > best) {
				best = aux;
				iBest = i;
			}				
		}
	}
	else {	
		for(int i = 1; i < this->getNumIndividuals(); i++) {
			aux = population.at(i).getFitness();
			if(aux < best) {
				best = aux;
				iBest = i;
			}				
		}	
	}
	
	return population.at(iBest);
	
}

void RealPopulation::copyPopulation(RealPopulation &b) {

	copyPop(population, b.population);

}

double RealPopulation::getSize() {

	return population.size();

}

const RealPopulation& RealPopulation::operator = (const RealPopulation & rhs) {
	
	if(this == &rhs)
		return *this;
		
	population.clear();
	
	this->size = rhs.size;
	this->numIndividuals = rhs.numIndividuals;
	for(int i = 0; i < numIndividuals; i++)
		this->population.push_back(rhs.population.at(i));		
		
	return *this;
	
}

void RealPopulation::setIndFitness(int index, double value) {

	population.at(index).setFitness(value);

}

void RealPopulation::evaluateAll(int n, double max, double min) {
	for(int i = 0; i < this->getNumIndividuals(); i++) {
			population[i].setFitness(double((random() / float(RAND_MAX)) * ((max-min)+1) + min));
		
	}
}

double RealPopulation::getFitnessSum() {
	
	double sum = 0.0;
	
	for(int i = 0; i < this->getNumIndividuals(); i++) {
		sum += population[i].getFitness();
	}
	
	return sum;
}

int RealPopulation::getWorstIndividualIndex(int optimizationType) {
    double worst, aux;
    int iWorst = 0;

    worst = population[0].getFitness();

    if(optimizationType == Operators::MIN) {
        for(int i = 0; i < this->getNumIndividuals(); i++) {
            aux = population.at(i).getFitness();
            if(aux > worst) {
                worst = aux;
                iWorst = i;
            }				
        }
    }
    else {	
        for(int i = 0; i < this->getNumIndividuals(); i++) {
            aux = population.at(i).getFitness();
            if(aux < worst) {
                worst = aux;
                iWorst = i;
            }				
        }	
    }

    return iWorst;
}

void RealPopulation::printPopulation() {

	for(int i = 0; i < this->getNumIndividuals(); i++)
		cout << population[i] << "\n";

}
