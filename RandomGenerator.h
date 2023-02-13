#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

class RandomGenerator
{
public:
	RandomGenerator(int idum);
	virtual ~RandomGenerator();
	void setSeed(int seed);
	double generate();
	

private:
	int idum;
	
};

#endif /*RANDOMGENERATOR_H_*/
