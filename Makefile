all: ag #model

Operators.o:
	g++ -I../../sundials/include/cvode -I../../sundials/include -I../../sundials/include/nvector -L ../../sundials/src/nvec_ser -L ../../sundials/src/cvode -c Operators.cpp -g3 -o Operators.o -lsundials_cvode -lsundials_nvecserial -lm -w
	#g++ -c Operators.cpp -g3 -o Operators.o

RealGenome.o:
	g++ -c RealGenome.cpp -g3 -o RealGenome.o

RealPopulation.o:
	g++ -c RealPopulation.cpp -g3 -o RealPopulation.o

main.o:
	../../mpich-1.2.7p1/bin/mpiCC -o main.o -g3 -c main_asinc.cpp -lm

calculos.o:
	g++ -w -c calculos.c -o calculos.o

ag: Operators.o RealGenome.o RealPopulation.o main.o calculos.o
	../../mpich-1.2.7p1/bin/mpiCC -static -fPIC -I../../sundials/include -I../../sundials/include/cvode -I../../sundials/include/nvector -L ../../sundials/src/cvode/.libs -L ../../sundials/lib -o ag -g3 Operators.o RealGenome.o RealPopulation.o main.o calculos.o -lm -lsundials_cvode -lsundials_nvecserial -L ../../sundials/src/nvec_ser/.libs

%o:%.cpp
	g++ -c -g3 $<


model: calculos.o
	g++ -I/opt/sundials/include -I/opt/sundials/include/cvode -I/opt/sundials/include/sundials -I/opt/sundials/include/nvector -L /opt/sundials/lib -w -O3 -o model teste.hpp  calculos.o -lm -lsundials_cvode -lsundials_nvecserial

	#Bondarenko_V5.2_SolveODE.hpp
