CC=gcc
F77=gfortran
LIBS=-lm -lgfortran

OPT1=-O2 -o
OPT2=-O2 -c

libplicroutines.a: plicroutines.o
	ar rv libplicroutines.a plicroutines.o

test: test.f90 mesh.o libplicroutines.a
	$(F77) $(OPT1) test_1 test.f90 mesh.o libplicroutines.a

all: test

plicroutines.o: plicroutines.f90
	$(F77) $(OPT2) plicroutines.f90

mesh.o: mesh.f90
	$(F77) $(OPT2) mesh.f90

clean:
	rm -f *.o
