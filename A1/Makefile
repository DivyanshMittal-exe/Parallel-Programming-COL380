CC = g++
CFLAGS = -std=c++17 -fopenmp -O2
ROOT_DIR:= /home/higgboson/Codes/Sem 6/Parallel-Programming-COL380/A1
LDFLAGS = -L. -lrary

all: main.o
	$(CC) $(CFLAGS) -o exec main.o $(LDFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

run: all
	./exec

comp:
	$(CC) $(CFLAGS) -o comp comp.cpp

clean:
	rm -f exec main.o comp

test_1: all
	./exec input output_test_1

test_2: all
	./exec input2 output_test_2