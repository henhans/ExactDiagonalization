all: main

LIBS = -L$LIBRARY_PATH -I$INCLUDE -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_intel_thread -lmkl_core -lm  -liomp5 #-llapack -lblas
CC = g++
CFLAGS = -Os -funroll-loops -Wall -pedantic

main: main.o matrix.o
	$(CC) $(CFLAGS) main.o matrix.o -o main $(LIBS)

main.o: main.cpp
	$(CC) -c main.cpp

matrix.o: matrix.cpp
	$(CC) -c matrix.cpp

clean: 
	rm -rf *.o main

