CC = g++
CFLAGS = -std=c++17 -O3

CXX = mpic++
CXXFLAGS = -O3

all:exec
	mpirun -np 4 ./exec

exec: main.o
	$(CXX) $(CXXFLAGS) -o exec main.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp

testcase_maker: testcase_maker.o
	$(CC) $(CFLAGS) -o testcase_maker testcase_maker.o

testcase_maker.o: testcase_maker.cpp
	$(CC) $(CFLAGS) -c -o testcase_maker.o testcase_maker.cpp

testcase: testcase_maker
	./testcase_maker sample

clean:
	rm -f main.o exec testcase_maker.o testcase_maker
