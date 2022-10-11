CC=g++
FLAGS=-std=c++17 -Wall -O2
LIBS=-lm -lgsl

main: main.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o stellar_wind.o utils.o
	$(CC) $(FLAGS) -o main main.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o stellar_wind.o utils.o $(LIBS)

main.o: gas_flow.hpp \
				godunov.hpp \
				riemann_solvers/hll.hpp \
				riemann_solvers/hllc.hpp \
				riemann_solvers/exact.hpp \
				stellar_wind.hpp \
				utils.hpp \
				main.cpp
	$(CC) $(FLAGS) -c main.cpp $(LIBS)

static_solution: static_solution.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o utils.o stellar_wind.o
	$(CC) $(FLAGS) -o static_solution static_solution.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o stellar_wind.o utils.o $(LIBS)

static_solution.o: gas_flow.hpp godunov.hpp riemann_solvers/hll.hpp riemann_solvers/hllc.hpp riemann_solvers/exact.hpp utils.hpp static_solution.cpp stellar_wind.hpp
	$(CC) $(FLAGS) -c static_solution.cpp $(LIBS)

supersonic_source: supersonic_source.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o utils.o
	$(CC) $(FLAGS) -o supersonic_source supersonic_source.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o utils.o $(LIBS)

supersonic_source.o: gas_flow.hpp godunov.hpp riemann_solvers/hll.hpp riemann_solvers/hllc.hpp riemann_solvers/exact.hpp utils.hpp supersonic_source.cpp
	$(CC) $(FLAGS) -c supersonic_source.cpp $(LIBS)

shu_osher: shu_osher.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o utils.o
	$(CC) $(FLAGS) -o shu_osher shu_osher.o godunov.o hll.o hllc.o exact.o gas_flow.o base.o utils.o $(LIBS)

shu_osher.o: gas_flow.hpp godunov.hpp riemann_solvers/hll.hpp riemann_solvers/hllc.hpp riemann_solvers/exact.hpp utils.hpp shu_osher.cpp
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

stellar_wind.o: stellar_wind.hpp stellar_wind.cpp gas_flow.hpp
	$(CC) $(FLAGS) -c stellar_wind.cpp $(LIBS)

utils.o: utils.hpp utils.cpp gas_flow.hpp
	$(CC) $(FLAGS) -c utils.cpp $(LIBS)

gas_flow.o: gas_flow.hpp gas_flow.cpp
	$(CC) $(FLAGS) -c gas_flow.cpp $(LIBS)

clean:
	rm -rf *.o main
