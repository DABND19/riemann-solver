CC=gcc
CFLAGS=-Wall -I.
LIBS=-lm -lgsl

all: main

main: main.o riemann_problem.o godunov.o
	$(CC) $(CFLAGS) -o main riemann_problem.o godunov.o main.o $(LIBS)

main.o: riemann_problem.h godunov.h main.c
	$(CC) $(CFLAGS) -c main.c $(LIBS)

godunov.o: riemann_problem.h godunov.h godunov.c
	$(CC) $(CFLAGS) -c godunov.c $(LIBS)

riemann_problem.o: riemann_problem.h riemann_problem.c
	$(CC) $(CFLAGS) -c riemann_problem.c $(LIBS)

clean:
	rm -fr *.o main
