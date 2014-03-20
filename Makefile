CC=gcc
CXX=g++

CFLAGS = -g -Wall
#CFLAGS = -O3 -Wall
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LIB = -lz -pthread

all: ngsDist


parse_args.o: parse_args.cpp ngsDist.hpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

shared.o: shared.cpp shared.hpp
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared.cpp

ngsDist: ngsDist.cpp parse_args.o shared.o
	$(CXX) $(CFLAGS) $(DFLAGS) ngsDist.cpp parse_args.o shared.o $(LIB) -o ngsDist

test:
	@cd examples/; sh ./test.sh 2> /dev/null; cd ../

clean:
	@rm -f *~ *.o ngsDist examples/testA_*
