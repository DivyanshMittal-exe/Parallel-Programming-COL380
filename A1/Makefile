CC = g++
CFLAGS = -std=c++17 -fopenmp -O2
LDFLAGS = -L. -lrary

all: main.o
	$(CC) $(CFLAGS) -o exec main.o $(LDFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

run: all
	./exec

clean:
	rm -f exec main.o
