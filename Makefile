CC = g++
CCFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -I InfoNest/cpp -I include -I celerite/cpp/include

default:
	make lib -C InfoNest/cpp
	$(CC) $(CCFLAGS) -c src/AR1.cpp
	$(CC) $(CCFLAGS) -c src/Model.cpp
	$(CC) $(CCFLAGS) -c src/Oscillation.cpp
	$(CC) $(CCFLAGS) -c src/WhiteNoise.cpp
	$(CC) $(CCFLAGS) -c src/main.cpp
	$(CC) -L InfoNest/cpp -o main *.o -larmadillo -linfonest
	rm -f *.o

