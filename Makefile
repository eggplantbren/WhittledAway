CC = g++
CCFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -I InfoNest/cpp -I include

default:
	make -C InfoNest/cpp
	$(CC) $(CCFLAGS) -c src/MyModel.cpp
	$(CC) $(CCFLAGS) -c src/main.cpp
	$(CC) -L InfoNest/cpp -o main *.o -linfonest
	rm -f *.o

