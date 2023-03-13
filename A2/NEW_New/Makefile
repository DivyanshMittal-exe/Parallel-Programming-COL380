CC = g++
CFLAGS = -std=c++17 -O3

CXX = mpic++
CXXFLAGS = -O3


all:exec
	mpirun -np 4 ./exec --taskid=1 --inputpath=test1.gra --headerpath=test1.dat --outputpath=sample_out --verbose=0 --startk=1 --endk=6 --p=0

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

test_base: exec
	mpirun -np 5 ./exec --taskid=1 --inputpath=0.gra --headerpath=0.dat --outputpath=o_me_0.txt --verbose=0 --startk=1 --endk=10 --p=0
#	mpirun -np 5 ./exec --taskid=1 --inputpath=./test_case/test0/test-input-0.gra --headerpath=./test_case/test0/test-header-0.gra --outputpath=./test_case/test0/o_me_0.txt --verbose=1 --startk=1 --endk=10 --p=0
	#diff ./test_case/test0/o_me_0.txt ./test_case/test0/output0.txt
#	mpirun -np 5 ./exec --taskid=1 --inputpath=test-input-1.gra --headerpath=test-header-1.dat --outputpath=o_me_1.txt --verbose=1 --startk=1 --endk=10 --p=0
#	diff test_case/test1/o_me_1.txt test_case/test1/output1.txt
#	mpirun -np 5 exec --taskid=1 --inputpath=test_case/test2/test-input-2.gra --headerpath=test_case/test2/test-header-2.dat --outputpath=test_case/test2/o_me_2.txt --verbose=1 --startk=1 --endk=8 --p=0
#	diff ./test_case/test2/o_me_2.txt ./test_case/test2/output2.txt
#	mpirun -np 5 ./exec --taskid=1 --inputpath=./test_case/test3/test-input-3.gra --headerpath=./test_case/test3/test-header-3.dat --outputpath=./test_case/test3/o_me_3.txt --verbose=1 --startk=1 --endk=5 --p=0
#	diff ./test_case/test3/o_me_3.txt ./test_case/test3/output3.txt
#	mpirun -np 5 ./exec --taskid=1 --inputpath=./test_case/test4/test-input-4.gra --headerpath=./test_case/test4/test-header-4.dat --outputpath=./test_case/test4/o_me_4.txt --verbose=1 --startk=1 --endk=6 --p=0
#	diff ./test_case/test4/o_me_4.txt ./test_case/test4/output4.txt
#	mpirun -np 5 ./exec --taskid=1 --inputpath=./test_case/test5/test-input-5.gra --headerpath=./test_case/test5/test-header-5.dat --outputpath=./test_case/test5/o_me_5.txt --verbose=1 --startk=1 --endk=8 --p=0
#	diff ./test_case/test5/o_me_5.txt ./test_case/test5/output5.txt