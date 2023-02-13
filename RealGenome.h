#ifndef REALGENOME_H_
#define REALGENOME_H_

#include <cstdio>
#include <iostream>

using namespace std;

class RealGenome
{
public:
	RealGenome(int size);
	RealGenome(){cromossome = NULL;};
	RealGenome(const RealGenome &copy);
	virtual ~RealGenome();	
	const int getSize() const;
	void setSize(int size);
	double *getCromossome() const;
	double getGen(int i) const;
	void setGen(int i, double val);
	void initialize(double lowerBound, double upperBound);
	void initialize(double *lowerBound, double *upperBound);
	void initializeZeros();
	double &operator [] (int i);
	double operator [] (int i) const;
	const double getFitness() const;
	void setFitness(double fit);	
	const RealGenome& operator = (const RealGenome & rhs);
	friend ostream &operator<<( ostream &, const RealGenome & );    
	
protected:
	int size;
	double *cromossome;
	void copyCromossome(double *dest, double *src);
	double fitness;
};

#endif /*REALGENOME_H_*/
