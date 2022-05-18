CC=g++
FLAGS=-std=c++17 -Wall -O2
LIBS=-lm -lgsl

shu_osher: shu_osher.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o
	$(CC) $(FLAGS) -o shu_osher shu_osher.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o $(LIBS)

shu_osher.o: gas_flow.hpp godunov.hpp riemann_solvers/hll.hpp riemann_solvers/hllc.hpp riemann_solvers/exact.hpp shu_osher.cpp
	$(CC) $(FLAGS) -c shu_osher.cpp $(LIBS)

godunov.o: gas_flow.hpp godunov.hpp godunov.cpp
	$(CC) $(FLAGS) -c godunov.cpp $(LIBS)

exact.o: gas_flow.hpp riemann_solvers/base.hpp riemann_solvers/exact.hpp riemann_solvers/exact.cpp
	$(CC) $(FLAGS) -c riemann_solvers/exact.cpp $(LIBS)

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
