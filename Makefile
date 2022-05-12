CC=g++
FLAGS=-std=c++17 -Wall -O2
LIBS=-lm -lgsl

all: main

main: main.o godunov.o hllc.o gas_flow.o base.o
	$(CC) $(FLAGS) -o main main.o godunov.o hllc.o gas_flow.o base.o $(LIBS)

main.o: gas_flow.hpp godunov.hpp riemann_solvers/hllc.hpp main.cpp
	$(CC) $(FLAGS) -c main.cpp $(LIBS)

godunov.o: gas_flow.hpp godunov.hpp godunov.cpp
	$(CC) $(FLAGS) -c godunov.cpp $(LIBS)

hllc.o: gas_flow.hpp riemann_solvers/base.hpp riemann_solvers/hllc.hpp riemann_solvers/hllc.cpp
	$(CC) $(FLAGS) -c riemann_solvers/hllc.cpp $(LIBS)

hll.o: gas_flow.hpp riemann_solvers/base.hpp riemann_solvers/hll.hpp riemann_solvers/hll.cpp
	$(CC) $(FLAGS) -c riemann_solvers/hll.cpp $(LIBS)

base.o: riemann_solvers/base.hpp riemann_solvers/base.cpp
	$(CC) $(FLAGS) -c riemann_solvers/base.cpp $(LIBS)

gas_flow.o: gas_flow.hpp gas_flow.cpp
	$(CC) $(FLAGS) -c gas_flow.cpp $(LIBS)

clean:
	rm -rf *.o main
