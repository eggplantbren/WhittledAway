CC = g++
CCFLAGS = -std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -Wall -Wextra -pedantic -I InfoNest/cpp -I include -I celerite/cpp/include -DWITH_LAPACK

default: infonest objects binaries

infonest:
	make lib -C InfoNest/cpp

objects: # Ones without a main() in them
	$(CC) $(CCFLAGS) -c src/AR1.cpp
	$(CC) $(CCFLAGS) -c src/Model.cpp
	$(CC) $(CCFLAGS) -c src/Oscillation.cpp
	$(CC) $(CCFLAGS) -c src/WhiteNoise.cpp

binaries:
	$(CC) $(CCFLAGS) -c src/main.cpp
	$(CC) -L InfoNest/cpp -o main *.o -linfonest
	rm -f main.o

	$(CC) $(CCFLAGS) -c test/test.cpp
	$(CC) -L InfoNest/cpp -o test/test *.o -linfonest
	rm -f test.o

	rm -f *.o

