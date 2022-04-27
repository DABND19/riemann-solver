CC=gcc
CFLAGS=-Wall -I. -O2
LIBS=-lm -lgsl

all: main

main: main.o riemann_problem.o godunov.o gas_flow.o
	$(CC) $(CFLAGS) -o main gas_flow.o riemann_problem.o godunov.o main.o $(LIBS)

main.o: gas_flow.h riemann_problem.h godunov.h main.c
	$(CC) $(CFLAGS) -c main.c $(LIBS)

godunov.o: gas_flow.h riemann_problem.h godunov.h godunov.c
	$(CC) $(CFLAGS) -c godunov.c $(LIBS)

riemann_problem.o: gas_flow.h riemann_problem.h riemann_problem.c
	$(CC) $(CFLAGS) -c riemann_problem.c $(LIBS)

gas_flow.o: gas_flow.h gas_flow.c
	$(CC) $(CFLAGS) -c gas_flow.c $(LIBS)

clean:
	rm -fr *.o main
