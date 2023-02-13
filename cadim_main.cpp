#include "RealPopulation.h"
#include "Operators.h"
#include <iostream>
#include <cmath>
#include <sstream>

#define SPACE " "
using namespace std;


#include"mpi.h"

#define IND_TAG 10
#define FIT_TAG 20

using namespace std;

MPI::Status status;	

void sendIndividual(RealGenome ind, int dest) {
   double *bits;
   int nBits = ind.getSize();
   bits = (double*)malloc(sizeof(double)*nBits);
   for(int i = 0; i < nBits; i++) {
            bits[i] = ind.getGen(i);
   }
   MPI::COMM_WORLD.Send(&bits[0], nBits, MPI::DOUBLE, dest, IND_TAG);
}

RealGenome receiveIndividual(int nBits) {	
   
   RealGenome b(nBits);
   double *bits;	
   bits = (double*)malloc(sizeof(double)*nBits);
   
   MPI::COMM_WORLD.Recv(&bits[0], nBits, MPI::DOUBLE, 0, IND_TAG, status);
            
   for(int i = 0; i < b.getSize(); i++) {
            b.setGen(i,bits[i]);		
   }	
   
   return b;
}

void sendFitness(double *fit, int nPac, int dest) {
 
   MPI::COMM_WORLD.Send(&fit[0], nPac, MPI::DOUBLE, dest, FIT_TAG);		
     
}

double *receiveFitness(int nPac, int src) {
   double *fit;	
   fit = (double*)malloc(sizeof(double)*nPac);
   MPI::COMM_WORLD.Recv(&fit[0], nPac, MPI::DOUBLE, src, FIT_TAG, status);
   return fit;
}

int main( int argc, char *argv[] ) {
      
   int myrank, nprocs;
   MPI::Init(argc, argv);    

   
   myrank = MPI::COMM_WORLD.Get_rank();
   nprocs = MPI::COMM_WORLD.Get_size(); 
      
   int n = 8;

   if(argc != 4) {
      fprintf(stderr,"To run the program use %s numGens nInd filename\n", argv[0]);
      exit(0);
   }

   int numGens = atoi(argv[1]);
   int nInd = atoi(argv[2]);
   char* filename = (char*) argv[3];
      
   int optType = Operators::MIN;
   int nSlaves = nprocs-1;
   int nLoad = nInd/nSlaves;
   
   double mutationRate = 0.25;
   double crossoverRate = 0.85;

   double *l = new double[n]; 
   double *r = new double[n]; 

         
   l[0] = 0.4067/1.5;    
   r[0] = 0.4067*1.5; 
   
   l[1] = 0.0;   
   r[1] = 1.5;
   
   l[2] = 0.16/1.5;   
   r[2] = 0.16*1.5;
   
   l[3] = 0.05/1.5;   
   r[3] = 0.05*1.5;
   
   l[4] = 0.078/1.5;   
   r[4] = 0.078*1.5;
   
   l[5] = 0.00575/1.5;   
   r[5] = 0.00575*1.5;
   
   l[6] = 13.0/1.5;   
   r[6] = 13.0*1.5; 
   
   l[7] = 0.1729/1.5; 
   r[7] = 0.1729*1.5; 
   
   double *fit;
               
   fit = (double*)malloc(nLoad*sizeof(double));	

   if ( myrank == 0 ) {
   
      float itime = MPI::Wtime();
      srand(time(NULL));	
      RealPopulation pop(nInd,n, l,r);
      
      RealPopulation newPop(nInd,n,l,r);
      RealGenome *parents;
      RealGenome *offspring;
      RealGenome best;
      double sum;
      int count = 0;	


      for(int slave = 1; slave <= nSlaves; slave++) {				
               for(int k = 0; k < nLoad; k++) {
                     sendIndividual(pop.getIndividual(count), slave);
                     count++;				
               } 				
      }						
      
      count = 0;
                     
      for(int slave = 1; slave <= nSlaves; slave++) {
      		printf("entrou\n");				
               fit = receiveFitness(nLoad,slave);						
               for(int k = 0; k < nLoad; k++) {											
                     pop.setIndFitness(count,fit[k]);
                     count++;				
               }			
      }
         
      Operators::savePopulation(pop,0);
      
      for(int gen = 1; gen < numGens; gen++) {       	
         
         Operators::quickSort(pop, 0, nInd-1, optType);
         best = pop.getIndividual(0);
         
         Operators::savePopulation(pop,gen);
         
         for(int j = 0; j < nInd; j+=2) {                
               parents = Operators::selectParentsByRank(pop);
               if(parents[0].getFitness() < parents[1].getFitness())
                     offspring = Operators::blendCrossover(parents[0],parents[1],crossoverRate);
               else
                     offspring = Operators::blendCrossover(parents[1],parents[0],crossoverRate);   

               Operators::mutation(offspring[0], mutationRate, l, r, gen, numGens);
               Operators::mutation(offspring[1], mutationRate, l, r, gen, numGens);                              

               newPop.setIndividual(j,offspring[0]);

               if((j+1) < nInd)
                     newPop.setIndividual(j+1,offspring[1]);
         }

         newPop.setIndividual(0,best);

         pop = newPop;

                     cout << gen << ": " << best << " - " << best.getFitness() <<  endl;			
                     
         count = 0;

         for(int i = 1; i <= nSlaves; i++) {				
            for(int k = 0; k < nLoad; k++) {
               sendIndividual(pop.getIndividual(count), i);
               count++;				
            } 				
         }

         count = 0;

         for(int slaves = 1; slaves <= nSlaves; slaves++) {				
            fit = receiveFitness(nLoad,slaves);			
            for(int k = 0; k < nLoad; k++) {									
               pop.setIndFitness(count,fit[k]);
               count++;				
            }			
         }			
      }

      Operators::quickSort(pop, 0, nInd-1, optType);
      Operators::savePopulation(pop,numGens);
      
      cout << numGens << ": "<< pop.getBestIndividual(optType) << " - "  << pop.getBestIndividual(optType).getFitness()  <<  endl;
               
      cout << "Time Elapsed (s): " << MPI::Wtime() - itime<< " s" << endl;
      cout << "Time Elapsed (h): " << (MPI::Wtime() - itime)/3600.0 << " h" << endl; 
   }
      
   else {
      int p = 0;
      RealGenome individual[nLoad];

      while(p < numGens) {
         for(int i = 0; i < nLoad; i++) {
            individual[i] = receiveIndividual(n);
            fit[i] = Operators::evaluateAllFrequencies(individual[i], filename);            
         }		
         sendFitness(fit,nLoad,0);
         p++;	
      }	
   }
   
   MPI::Finalize();
}
