CC = g++
CFLAGS = -std=c++17 -O3

CXX = mpic++
CXXFLAGS = -O3


all:exec
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test3/test-input-3.gra --headerpath=test_case/test3/test-header-3.dat --outputpath=test_case/test3/o_me.txt --verbose=0 --startk=1 --endk=5 --p=0

exec: main.o
	$(CXX) $(CXXFLAGS) -o exec main.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o main.cpp

testcase_maker: testcase_maker.o
	$(CC) $(CFLAGS) -o testcase_maker testcase_maker.o

testcase_maker.o: testcase_maker.cpp
	$(CC) $(CFLAGS) -c -o testcase_maker.o testcase_maker.cpp

testcase: testcase_maker
	./testcase_maker sample3

clean:
	rm -f main.o exec testcase_maker.o testcase_maker

test_css: exec
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test0/test-input-0.gra --headerpath=test_case/test0/test-header-0.dat --outputpath=test_case/test0/o_me.txt --verbose=1 --startk=1 --endk=10 --p=0
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test1/test-input-1.gra --headerpath=test_case/test1/test-header-1.dat --outputpath=test_case/test1/o_me.txt --verbose=1 --startk=1 --endk=10 --p=0
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test2/test-input-2.gra --headerpath=test_case/test2/test-header-2.dat --outputpath=test_case/test2/o_me.txt --verbose=1 --startk=1 --endk=8 --p=0
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test3/test-input-3.gra --headerpath=test_case/test3/test-header-3.dat --outputpath=test_case/test3/o_me.txt --verbose=1 --startk=1 --endk=5 --p=0
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test4/test-input-4.gra --headerpath=test_case/test4/test-header-4.dat --outputpath=test_case/test4/o_me.txt --verbose=1 --startk=1 --endk=6 --p=0
	mpirun -np 8 ./exec --taskid=1 --inputpath=test_case/test5/test-input-5.gra --headerpath=test_case/test5/test-header-5.dat --outputpath=test_case/test5/o_me.txt --verbose=1 --startk=1 --endk=8 --p=0

test_base: exec
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test0/test-input-0.gra --headerpath=test_case/test0/test-header-0.dat --outputpath=test_case/test0/o_me.txt --verbose=1 --startk=1 --endk=10 --p=0
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test1/test-input-1.gra --headerpath=test_case/test1/test-header-1.dat --outputpath=test_case/test1/o_me.txt --verbose=1 --startk=1 --endk=10 --p=0
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test2/test-input-2.gra --headerpath=test_case/test2/test-header-2.dat --outputpath=test_case/test2/o_me.txt --verbose=1 --startk=1 --endk=8 --p=0
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test3/test-input-3.gra --headerpath=test_case/test3/test-header-3.dat --outputpath=test_case/test3/o_me.txt --verbose=1 --startk=1 --endk=5 --p=0
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test4/test-input-4.gra --headerpath=test_case/test4/test-header-4.dat --outputpath=test_case/test4/o_me.txt --verbose=1 --startk=1 --endk=6 --p=0
	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test5/test-input-5.gra --headerpath=test_case/test5/test-header-5.dat --outputpath=test_case/test5/o_me.txt --verbose=1 --startk=1 --endk=8 --p=0
#	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test6/test-input-6.gra --headerpath=test_case/test6/test-header-6.dat --outputpath=test_case/test6/o_me.txt --verbose=0 --startk=1 --endk=10 --p=0
#	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test7/test-input-7.gra --headerpath=test_case/test7/test-header-7.dat --outputpath=test_case/test7/o_me.txt --verbose=0 --startk=1 --endk=10 --p=0
	#mpirun -np 6 ./exec --taskid=1 --inputpath=test_case/test8/test-input-8.gra --headerpath=test_case/test8/test-header-8.dat --outputpath=test_case/test8/o_me.txt --verbose=1 --startk=12 --endk=29 --p=0
#	mpirun -np 3 ./exec --taskid=1 --inputpath=test_case/test9/test-input-9.gra --headerpath=test_case/test9/test-header-9.dat --outputpath=test_case/test9/o_me.txt --verbose=0 --startk=1 --endk=10 --p=0
	#diff test_case/test0/o_me.txt test_case/test0/output0_verbose.txt
	#diff test_case/test1/o_me.txt test_case/test1/output1_verbose.txt
	#diff test_case/test2/o_me.txt test_case/test2/output2_verbose.txt
#	diff test_case/test3/o_me.txt test_case/test3/output3_verbose.txt
#	diff test_case/test4/o_me.txt test_case/test4/output4_verbose.txt
#	diff test_case/test5/o_me.txt test_case/test5/output5_verbose.txt
