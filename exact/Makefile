CC = gcc
CPPC = g++
CFLAGS = -I. -O3 -Wall -Wextra -std=c99 -lgmp -lm
CXXFLAGS = -I. -O3 -Wall -Wextra -std=c++17 -lgmp -lm

all: test.exe bench.exe

test.exe : test.o libexactWigner.a
	$(CC) -o test.exe test.o -L. -lexactWigner $(CFLAGS)

bench.exe : bench.o libexactWigner.a
	$(CPPC) -o bench.exe bench.o -L. -lexactWigner $(CXXFLAGS)

test.o : test.c
	$(CC) -c test.c $(CFLAGS)

bench.o : bench.cpp
	$(CPPC) -c bench.cpp $(CXXFLAGS)

libexactWigner.a : exactWigner.o
	ar rcs libexactWigner.a exactWigner.o
	ranlib libexactWigner.a

exactWigner.o : exactWigner.c
	$(CC) -c exactWigner.c $(CFLAGS)

clean:
	rm -f *.o *.a *.exe
	rm -f test.exe
	rm -f libexactWigner.a
	rm -f exactWigner.o
	rm -f test.o