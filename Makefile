CC=g++
FLAGS=-std=c++17 -Wall -O2
LIBS=-lm -lgsl

main: main.o \
			hll.o \
			gas_flow.o \
			base_godunov.o \
			base.o
	$(CC) $(FLAGS) -o main main.o hll.o gas_flow.o base_godunov.o base.o $(LIBS)

main.o: main.cpp \
				gas_flow.hpp \
				base_godunov.hpp \
				riemann_solvers/hll.hpp \
				shared.hpp
	$(CC) $(FLAGS) -c main.cpp $(LIBS)

base_godunov.o: base_godunov.hpp \
								base_godunov.cpp \
								gas_flow.hpp \
								riemann_solvers/base.hpp
	$(CC) $(FLAGS) -c base_godunov.cpp $(LIBS)

hll.o: gas_flow.hpp \
			 riemann_solvers/base.hpp \
			 riemann_solvers/hll.hpp \
			 riemann_solvers/hll.cpp
	$(CC) $(FLAGS) -c riemann_solvers/hll.cpp $(LIBS)

base.o: riemann_solvers/base.hpp riemann_solvers/base.cpp
	$(CC) $(FLAGS) -c riemann_solvers/base.cpp $(LIBS)

gas_flow.o: gas_flow.hpp \
						gas_flow.cpp
	$(CC) $(FLAGS) -c gas_flow.cpp $(LIBS)

clean:
	rm -rf *.o main
