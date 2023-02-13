#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;


#include"mpi.h"

#define IND_TAG 10
#define FIT_TAG 20

using namespace std;

MPI::Status status;	




int main( int argc, char *argv[] ) {
      
   int myrank, nprocs;
   MPI::Init(argc, argv);    

   
   myrank = MPI::COMM_WORLD.Get_rank();
   nprocs = MPI::COMM_WORLD.Get_size(); 
      
   int n = 8;

}
