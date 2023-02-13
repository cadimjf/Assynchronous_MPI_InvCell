#include "Operators.h"
#include <cmath>
#include <cstdio>
#include <sstream>
#include<vector>
#include <fstream>
#include <stdexcept>

#include "calculos.h"
#include "Bondarenko_V5.2_SolveODE.hpp"

using namespace std;

int Operators::randI( int min = 0 , int max = 1 ) {
	return int((random() / float(RAND_MAX)) * ((max-min)+1) + min);
}

double Operators::randD( double min = 0.0, double max = 1.0 ) {
	return random() / double(RAND_MAX) * (max-min) + min;  
}

/*double Operators::randD( double min = 0.0, double max = 1.0 ) {
	return ((double)(rand()%1000)/1000.0)*(max-min) + min;  
}*/

RealPopulation Operators::restorePopulation(string fileName, int nInd, int nP) {

	RealPopulation pop(nInd,nP);	
	FILE *f = fopen(fileName.c_str(), "r");
	double aux;

	if(!f) {
		cerr << "Fail to open population File " << fileName << "\n" << "Returning zero population!\n";
		return pop;
	}

	for(int i = 0; i < nInd; i++) {
		for(int j = 0; j < nP; j++) {
			fscanf(f,"%lf",&aux);
			pop.getIndividual(i).setGen(j,aux);
		}
	}

	return pop;

}

int Operators::coin() {

	double r = randD();

	if(r < 0.5) {	
		return 0;	
	}	
	else return 1;

}

inline void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

inline int w(double v, double v_) {

	double d = v - v_;

	if( (d < 0.0) || (fabs(v + 60.0) < 0.001) )
		return 1000;
	else
		return 1;

}
inline vector<double> readData(string fileName, double t) {

   vector<string> tokens;
   fstream file_op(fileName.c_str(),ios::in);
   string str;
   vector<double> data;
   double time;
   
   if(!file_op.is_open()) {
      cout << "Error opening file" << endl;
      exit(1);
   }
   
   while(getline(file_op,str)) {
      tokenize(str,tokens," ");
      time = atof(tokens[0].c_str());        
      if(time >= t) {			     
            data.push_back(atof(tokens[1].c_str()));	        	        
      }
      tokens.clear();                            
   }
   
   return data;
}

//inline vector<double> readData(string fileName) {

//	vector<string> tokens;
//	fstream file_op(fileName.c_str(),ios::in);
//	string str;
//	vector<double> data;
//	double time;

//    if(!file_op.is_open()) {
//        cout << "Error opening file" << endl;
//        exit(1);
//    }

//	while(getline(file_op,str)) {
//		tokenize(str,tokens," ");
//		time = atof(tokens[0].c_str());        
//		data.push_back(atof(tokens[1].c_str()));	        	        
//		tokens.clear();                            
//	}

//	return data;
//}

inline vector<vector<double> > readData(string configFile) {

	fstream file_op(configFile.c_str(),ios::in);
	
	string str;
	string str1;	
	vector<vector<double> > data;
    

    if(!file_op.is_open()) {
        cout << "Error opening file " << configFile << endl;
        exit(1);
    }

	while(getline(file_op,str)) {
	
        vector <double> aux;	
		fstream data_file(str.c_str(),ios::in);
	
	    if(!data_file.is_open()) {
            cout << "Error opening file " << str << endl;
            exit(1);
        }
		
		while(getline(data_file,str1)) {
		
		    aux.push_back(atof(str1.c_str()));    
		
		}
		
		data_file.close();
		data.push_back(aux);
		
	}

	return data;
}

inline vector<double> readExpData(string fileName) {

    fstream file_op(fileName.c_str(),ios::in);
    string str;
    vector<double> data;

    if(file_op.is_open()) { 
        while(getline(file_op,str)) {                        
            data.push_back(atof(str.c_str()));                                           
        }
    }
    else {
        cout << "Failed opening " << fileName << endl;
        exit(0);
    }

    return data;
}

// inline vector<double> readExpData(string fileName) {
// 
// 	vector<string> tokens;
// 	fstream file_op(fileName.c_str(),ios::in);
// 	string str;
// 	vector<double> data;
// 
// 	if(file_op.is_open()) {	
// 		while(getline(file_op,str)) {
// 			tokenize(str,tokens," ");        
// 			data.push_back(atof(tokens[1].c_str()));	        	        
// 			tokens.clear();                            
// 		}
// 	}
// 	else {
// 		cout << "Failed opening " << fileName << endl;
// 	}
// 
// 	return data;
// }

inline double errorFunction(double *cond, vector<double> dd) { 

   Solveode *EDO = new Solveode();

   FILE *f1;

   double *t,
          *vm,
          *cai,
          *nai,
          *I_Ktof,
          *I_Kur,
          *I_Kss,
          *I_Ks, 
          *I_Kr, 
          *I_CaL;

   double IK_total[35000];

   double g_Kto_f,
          g_Kto_s,
          g_Kur,
          g_Kss,
          g_CaL,
          g_Na,
          g_Ks,
          g_Kr,
          dt;

   // medidas
   // PA
   double repouso,
          pico,
          amplitude,
          max_dVdt,
          APD30,
          APD50,
          APD90;

   double repouso_exp,
          pico_exp,
          amplitude_exp,
          max_dVdt_exp,
          APD30_exp,
          APD50_exp,
          APD90_exp;

   double maxIK[12],
          minIK[12];

   double maxIK_exp[12],
          minIK_exp[12];

   double maxICa[10],
          minICa[10];

   double maxICa_exp[10],
          minICa_exp[10];

   int i,
       j,
       n,
       v_step;

   char arq[10],
        vstep[2];

   double error;

   // Parâmetros
   g_Kto_f = cond[0];
   g_Kto_s = cond[1];
   g_Kur   = cond[2];
   g_Kss   = cond[3];
   g_Kr    = cond[4];
   g_Ks    = cond[5];
   g_Na    = cond[6];
   g_CaL   = cond[7];

   // Discretização
   dt = 0.1;

/*--------------------------------Resolvendo o modelo--------------------------------*/
   EDO->setFreeVariable(dt);

   // Pacing
   n  = 100000;

   // Setando os parâmetros
   EDO->setParameters(54, g_Kto_f);
   EDO->setParameters(55, g_Kto_s);
   EDO->setParameters(57, g_Kur);
   EDO->setParameters(58, g_Kss);
   EDO->setParameters(59, g_Kr);
   EDO->setParameters(56, g_Ks);
   EDO->setParameters(51, g_Na);
   EDO->setParameters(34, g_CaL);

   EDO->setParameters(74, 0); // protocolo de pacing
   EDO->solveCVODE(1, n);

   t  = EDO->getIndependentVar();
   vm = EDO->getSolution(0);
   cai = EDO->getSolution(1);
   nai = EDO->getSolution(18);

   // Medidas
   amplitude = Amplitude(10000, &vm[90000], &pico, &repouso);
   max_dVdt = Max_dVdt(10000, dt, &vm[90000]);
   APD(10000, dt, &vm[90000], &t[90000], &APD30, 30);
   APD(10000, dt, &vm[90000], &t[90000], &APD50, 50);
   APD(10000, dt, &vm[90000], &t[90000], &APD90, 90);


   // Protocolo de voltage step da IK
   v_step = -40;
   for(i = 0; i < 12; i++) {
      n  = 35000;

      EDO->setParameters(0, 0.0); //reinicia o tempo com 0
      EDO->reInitCVODE(); // reinicia o CVode
      EDO->setVariablesFromFile("condicoes_iniciais");

      // Setando os parâmetros
      EDO->setParameters(54, g_Kto_f);
      EDO->setParameters(55, g_Kto_s);
      EDO->setParameters(57, g_Kur);
      EDO->setParameters(58, g_Kss);
      EDO->setParameters(59, g_Kr);
      EDO->setParameters(56, g_Ks);
      EDO->setParameters(51, g_Na);
      EDO->setParameters(34, g_CaL);

      EDO->setParameters(73, v_step);
      EDO->setParameters(74, 1);   // protocolo da IK
      EDO->solveCVODE(1, n);

      t  = EDO->getIndependentVar();
   
      I_Ktof = EDO->getSolution(100);  
      I_Kur  = EDO->getSolution(101);  
      I_Kss  = EDO->getSolution(102);  
      I_Ks   = EDO->getSolution(103);  
      I_Kr   = EDO->getSolution(104);  

      for(j = 0; j < n; j++)
         IK_total[j] = I_Ktof[j] + I_Kur[j] + I_Kss[j] + I_Ks[j] + I_Kr[j]; 

      v_step += 10;

      // Medidas
      ind_max(n, IK_total, &maxIK[i]);
      ind_min(n, IK_total, &minIK[i]);
   }



   // Protocolo de voltage step da ICa
   v_step = -40;
   for(i = 0; i < 10; i++) {
      n  = 4000;

      EDO->setParameters(0, 0.0);
      EDO->reInitCVODE();
      EDO->setVariablesFromFile("condicoes_iniciais");

      // Setando os parâmetros
      EDO->setParameters(54, g_Kto_f);
      EDO->setParameters(55, g_Kto_s);
      EDO->setParameters(57, g_Kur);
      EDO->setParameters(58, g_Kss);
      EDO->setParameters(59, g_Kr);
      EDO->setParameters(56, g_Ks);
      EDO->setParameters(51, g_Na);
      EDO->setParameters(34, g_CaL);

      EDO->setParameters(73, v_step);
      EDO->setParameters(74, 2);   // protocolo da ICa
      EDO->solveCVODE(1, n);

      t = EDO->getIndependentVar();
   
      I_CaL = EDO->getSolution(106);

      v_step += 10;

      // Medidas
      ind_max(n, I_CaL, &maxICa[i]);
      ind_min(n, I_CaL, &minICa[i]);
   }

/*-----------------------------------------------------------------------------------*/

   APD30_exp     = dd[0];
   APD50_exp     = dd[1];
   APD90_exp     = dd[2];
   max_dVdt_exp  = dd[3];
   pico_exp      = dd[4];
   repouso_exp   = dd[5];
   amplitude_exp = dd[6];

   for(int i = 7; i < 19; i++)
      maxIK_exp[i-7] = dd[i];

   for(int i = 19; i < 31; i++)
      minIK_exp[i-19] = dd[i];

   for(int i = 31; i < 41; i++)
      maxICa_exp[i-31] = dd[i];

   for(int i = 41; i < 51; i++)
      minICa_exp[i-41] = dd[i];

   error = 0.0;

   error += fabs(repouso_exp - repouso)/fabs(repouso_exp);
   error += fabs(amplitude_exp - amplitude)/fabs(amplitude_exp);
   error += fabs(max_dVdt_exp - max_dVdt)/fabs(max_dVdt_exp);
   error += fabs(APD30_exp - APD30)/fabs(APD30_exp);
   error += fabs(APD50_exp - APD50)/fabs(APD50_exp);
   error += fabs(APD90_exp - APD90)/fabs(APD90_exp);

   for(int i = 0; i < 12; i++) {
      error += fabs(maxIK_exp[i] - maxIK[i])/fabs(maxIK_exp[i]);
      error += fabs(minIK_exp[i] - minIK[i])/fabs(minIK_exp[i]);
   }

   for(int i = 0; i < 10; i++) {
      error += fabs(maxICa_exp[i] - maxICa[i])/fabs(maxICa_exp[i]);
      error += fabs(minICa_exp[i] - minICa[i])/fabs(minICa_exp[i]);
   }

   delete EDO;

   return error;
} 

double Operators::evaluateAllFrequencies(RealGenome ind, char* filename) {
   vector<double> desirableData;
   desirableData = readExpData(filename);

   return errorFunction(ind.getCromossome(), desirableData);
}

inline double getDistance(double *v1, double *v2, int n) {

	double distance = 0.0;

	for(int i = 0; i < n; i++) {	
		distance += fabs(v1[i] - v2[i]);	
	}

	return distance;
}

RealGenome * Operators::selectParents(RealPopulation &b, int optmizationType) {

	int rand = randI(0, b.getNumIndividuals()-1);
	RealGenome ind1, ind2;
	RealGenome *parents = new RealGenome[2];

	ind1 = b.getIndividual(rand);

	rand = randI(0, b.getNumIndividuals()-1);

	ind2 = b.getIndividual(rand);

	if(optmizationType == MAX) {

		if(ind1.getFitness() >= ind2.getFitness())
			parents[0] = ind1;
		else
			parents[0] = ind2;

	} else {

		if(ind1.getFitness() <= ind2.getFitness())
			parents[0] = ind1;
		else
			parents[0] = ind2;

	}

	rand = randI(0, b.getNumIndividuals()-1);		

	ind1 = b.getIndividual(rand);

	rand = randI(0, b.getNumIndividuals()-1);

	ind2 = b.getIndividual(rand);        

	if(optmizationType == MAX) {

		if(ind1.getFitness() >= ind2.getFitness())
			parents[1] = ind1;
		else
			parents[1] = ind2;

	} else {

		if(ind1.getFitness() <= ind2.getFitness())
			parents[1] = ind1;
		else
			parents[1] = ind2;

	}

	return parents;	
}

RealGenome *Operators::selectParentsByRoullette(RealPopulation &pop, double sum) {

	double tempSum = 0.0;
	double value = randD()*(1.0/sum);
	int index = 0;
	int popSize = pop.getNumIndividuals();
	RealGenome *parents = new RealGenome[2];

	while ( (tempSum < value) && (index < (popSize-1) )) {
		tempSum += pop.getIndividual(index).getFitness();
		index++;
	}

	parents[0] = pop.getIndividual(index);

	tempSum = 0.0;
	value = randD()*(1.0/sum);
	index = 0;

	while ( (tempSum < value) && (index < (popSize-1) )) {
		tempSum += pop.getIndividual(index).getFitness();
		index++;
	}

	parents[1] = pop.getIndividual(index);

	return parents;
}

RealGenome * Operators::selectParentsByRank(RealPopulation &pop) {

	double rand;
	int indPop;
	RealGenome *parents = new RealGenome[2];

	indPop = (int)(pop.getNumIndividuals() * (1.5 - sqrt(2.25 - 2 * random())));

	if (indPop < 0) {
		indPop = 0;             
	}

	parents[0] = pop.getIndividual(indPop);

	indPop = (int)(pop.getNumIndividuals() * (1.5 - sqrt(2.25 - 2 * random())));

	if (indPop < 0) {
		indPop = 0;             
	}

	parents[1] = pop.getIndividual(indPop);

	return parents;

}

RealGenome Operators::crossover(RealGenome *parents, double probCrossover) {

	int j;
	double r = randD();
	int n = parents[0].getSize();
	RealGenome offspring(n);


	if(r < probCrossover) {		
		for(j = 0; j < n; j++) {
			if(parents[0].getGen(j) != parents[1].getGen(j)) {
				offspring.setGen(j,(parents[0].getGen(j) + parents[1].getGen(j))/2.0); 
			}
			else
				offspring.setGen(j,parents[0].getGen(j));
		}

		return offspring;		
	} 
	else
		return parents[0];
}

inline double generateBeta(double alfa) {

	double g = (double)((unsigned int)random()/(double)(RAND_MAX));
	double e = (int)random()%2;
	e = (float)e + g - alfa;			
	return e;
}

RealGenome * Operators::blendCrossover(RealGenome x, RealGenome y, double probCrossover) {

	RealGenome * offspring = new RealGenome[2];
	double prob = randD();
	int n = x.getSize();

	double a = 0.6,b = 0.7,di;
	double u;

	if (prob < probCrossover) {

		offspring[0].setSize(n);
		offspring[1].setSize(n);

		offspring[0].initializeZeros();
		offspring[1].initializeZeros();		

		for(int i = 0; i < n; i++) {

			di = fabs(x.getGen(i)-y.getGen(i));

			if (x.getGen(i) <= y.getGen(i)) {

				u = randD(x.getGen(i)-a*di, y.getGen(i)+b*di);
				offspring[0].setGen(i,u);

				u = randD(x.getGen(i)-a*di, y.getGen(i)+b*di);
				offspring[1].setGen(i,u);
			}
			else {
				u = randD(y.getGen(i)-b*di, x.getGen(i)+a*di);
				offspring[0].setGen(i,u);

				u = randD(y.getGen(i)-b*di, x.getGen(i)+a*di);
				offspring[1].setGen(i,u);
			}
		}

		return offspring;
	}
	else {
		offspring[0] = x;
		offspring[1] = y;
		return offspring;
	}
}

inline double delta(double t, double y, int T) {

	double b = 5.0;

	return y*(1.0 - pow(Operators::randD(),pow(1.0-(t/T),b)));

}

void Operators::mutation(RealGenome &ind, double probMutation, double *LB, double *UB, int t, int T)
{
	double ran;
	int i, d;
	int l = ind.getSize();	

	for(i = 0; i < l; i++){
		double r = randD();
		if(r < probMutation) {
			if(randI() == 0)
				ind.setGen(i,ind.getGen(i) + delta(t,UB[i]-ind.getGen(i),T));
			else
				ind.setGen(i,ind.getGen(i) - delta(t,ind.getGen(i)-LB[i],T));
		}
	}
}	

void Operators::mutation(RealGenome &ind, double probMutation, double LB, double UB, int t, int T)
{
	double ran;
	int i, d;
	int l = ind.getSize();	

	for(i = 0; i < l; i++){
		double r = randD();
		if(r < probMutation) {
			if(randI() == 0)
				ind.setGen(i,ind.getGen(i) + delta(t,UB-ind.getGen(i),T));
			else
				ind.setGen(i,ind.getGen(i) - delta(t,ind.getGen(i)-LB,T));
		}
	}
}	

void Operators::quickSort(RealPopulation &pop, int lowerBound, int supBound, int optType ) {

	int l,r;
	RealGenome pivot,aux;
	l = lowerBound;
	r = supBound;

	pivot = pop.getIndividual((int)((lowerBound+supBound)/2));    

	while(l <= r) 
	{
		if (optType == MAX) {

			while (pop.getIndividual(l).getFitness() > pivot.getFitness()) 
				l = l + 1;

			while (pop.getIndividual(r).getFitness() < pivot.getFitness()) 
				r = r - 1;
		}
		else
		{
			while (pop.getIndividual(l).getFitness() < pivot.getFitness()) 
				l = l + 1;

			while (pop.getIndividual(r).getFitness() > pivot.getFitness()) 
				r = r - 1;
		}

		if (l <= r) {

			aux = pop.getIndividual(l);
			pop.setIndividual(l,pop.getIndividual(r));
			pop.setIndividual(r,aux);

			l = l + 1;
			r = r - 1;
		}
	}

	if ((r-lowerBound) > 0) 
		quickSort(pop,lowerBound,r, optType);
	if ((supBound - l) > 0) 
		quickSort(pop,l,supBound, optType);

}

void Operators::saveFitnessMap(RealPopulation &pop, double *bestParams, int generation, int n, double *l, double *r ) {

	stringstream fileName(stringstream::in | stringstream::out);
	double dist;

	fileName << "gens/Generation " << generation;

	fstream file_op;

	file_op.open (fileName.str().c_str(), fstream::out);

	for(int i = 0; i < pop.getNumIndividuals(); i++) {
		dist = getDistance(pop.getIndividual(i).getCromossome(),bestParams,n);
		file_op << pop.getIndividual(i).getFitness() << " " << dist << "\n";
	}	

	file_op.close();
}

void Operators::savePopulation(RealPopulation &pop, int gen) {

	stringstream fileName1(stringstream::in | stringstream::out);
        stringstream fileName2(stringstream::in | stringstream::out);

	fileName1 << "population.txt";
        fileName2 << "fitness.txt";

	fstream file_op1;
        fstream file_op2;

	file_op1.open (fileName1.str().c_str(), fstream::out | fstream::app);
        file_op2.open (fileName2.str().c_str(), fstream::out | fstream::app);

	file_op1 << "Generation " << gen << ": " << "\n";

	for(int i = 0; i < pop.getNumIndividuals(); i++) {
		file_op1 << pop.getIndividual(i) << " - " << pop.getIndividual(i).getFitness() << "\n";
                file_op2 << gen << pop.getIndividual(i).getFitness() << "\n";
	}	

	file_op1 << "------------------------------------------------------\n";

	file_op1.close();
        file_op2.close();
}
