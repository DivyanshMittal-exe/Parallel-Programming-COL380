CC = g++
CFLAGS = -std=c++17 -fopenmp -O3
ROOT_DIR:= /home/higgboson/Codes/Sem 6/Parallel-Programming-COL380/A1
LDFLAGS = -L. -lrary -Wl,-rpath=.


all: main.o
	export LD_LIBRARY_PATH=./:$LB_LIBRARY_PATH
	$(CC) $(CFLAGS) -o exec main.o $(LDFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

run: all
	./exec

comp: comp.o
	$(CC) $(CFLAGS) -o comp comp.o

comp.o:comp.cpp
	$(CC) $(CFLAGS) -c comp.cpp


clean:
	rm -f exec main.o comp comp.o main main.o test_case_gen test_case_gen.o main_test.o exec_test

test_1: all
	./exec input output_test_1

test_2: all
	./exec input2 output_test_2


test_case_gen: test_case_gen.o
	$(CC) $(CFLAGS) -o test_case_gen test_case_gen.o $(LDFLAGS)

test_case_gen.o:test_case_gen.cpp
	$(CC) $(CFLAGS) -c test_case_gen.cpp