#include "RealGenome.h"
#include "RandomGenerator.h"
#include "Operators.h"

RealGenome::RealGenome(int size) {

	this->size = size;
	cromossome = new double[size];

}

RealGenome::RealGenome(const RealGenome &copy) {	
	this->size = copy.size;
	this->fitness = copy.fitness;
	this->cromossome = new double[size];
	copyCromossome(cromossome, copy.cromossome);
}

RealGenome::~RealGenome() {
	delete[] cromossome;
}

double * RealGenome::getCromossome() const {
	return this->cromossome;
}

double &RealGenome::operator [] (int i)  {

	if(i < this->size)
		return cromossome[i];
	else {
		cout << "Genome index out of bounds: " << i << endl;
		exit(0);
	}

}

double RealGenome::operator [] (int i) const {

	if(i < this->size)
		return cromossome[i];
	else {
		cout << "Genome index out of bounds: " << i << endl;
		exit(0);
	}

}

const double RealGenome::getFitness() const {

	return this->fitness;

}

void RealGenome::setFitness(double fit) {

	this->fitness = fit;

}

void RealGenome::copyCromossome(double *dest,double *src) {
	
	for(int i = 0; i < size; i++) {
	
		dest[i] = src[i];
	
	}

}

const RealGenome& RealGenome::operator = (const RealGenome & rhs){
	
	if(this == &rhs)
		return *this;
	
	if(cromossome != NULL) {
		delete[] cromossome;
	}
	
	this->size = rhs.size;	
	this->fitness = rhs.fitness;
	cromossome = new double[size];
	copyCromossome(cromossome, rhs.cromossome);	
	
	return *this;
}

const int RealGenome::getSize() const {

	return this->size;

}

void RealGenome::setSize(int size) {

	this->size = size;

}

ostream &operator<<( ostream &output, const RealGenome &a ) {

    for ( int i = 0; i < a.size; i++ ) {
    	output << a.cromossome[ i ] << " ";

     }

     return output;
}

double RealGenome::getGen(int i) const {

	return cromossome[i];

}

void RealGenome::setGen(int i, double val) {

	cromossome[i] = val;

}

void RealGenome::initialize(double lowerBound, double upperBound) {

	if(cromossome == NULL)
		cromossome = new double[size];

	for(int i = 0; i < size; i++)
		cromossome[i] = Operators::randD(lowerBound, upperBound); 

}

void RealGenome::initialize(double *lowerBound, double *upperBound) {

	if(cromossome == NULL)
		cromossome = new double[size];

	for(int i = 0; i < size; i++)
		cromossome[i] = Operators::randD(lowerBound[i], upperBound[i]); 

}

void RealGenome::initializeZeros() {

	if(cromossome == NULL)
		cromossome = new double[size];

	for(int i = 0; i < size; i++)
		cromossome[i] = 0.0; 

}

