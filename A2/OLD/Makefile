CC = g++
CFLAGS = -std=c++17 -O3

CXX = mpic++
CXXFLAGS = -Wall -Wextra -O3

all:exec
	mpirun -np 3 ./exec ../test_case/test0/test-input-0.gra ../test_case/test0/test-header-0.dat

exec: main.o
	$(CXX) $(CXXFLAGS) -o exec main.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp

testcase_maker: testcase_maker.o
	$(CC) $(CFLAGS) -o testcase_maker testcase_maker.o

testcase_maker.o: testcase_maker.cpp
	$(CC) $(CFLAGS) -c -o testcase_maker.o testcase_maker.cpp

test_case: testcase_maker
	./testcase_maker

clean:
	rm -f main.o exec test_casemaker.o testcase_maker
