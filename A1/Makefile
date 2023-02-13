CC = g++
CFLAGS = -std=c++17 -g -fopenmp -O3
ROOT_DIR:= /home/higgboson/Codes/Sem 6/Parallel-Programming-COL380/A1
LDFLAGS = -L. -lrary

all: main.o
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




analyse: test_case_gen
#	$(CC) $(CFLAGS) -c -Dtime_test main.cpp -o main_test.o
#	$(CC) $(CFLAGS) -o exec_test main_test.o $(LDFLAGS)
	#exec 3> results.txt \
	for i in {10..18}; do \
	  ./test_case_gen $i; \
	  for j in {1,2,4,8,16}; do \
		real_time=$(echo $(time OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm) | awk '{print $2}') \
#		rm output_mmmm; \
		echo $i $j $real_time; \
	  done; \
	  rm test_case_$i; \
	done

	# Close the file
	#exec 3>&-


	#perf	stat ./exec_test test_case_10 output_10
	#perf	stat ./exec_test test_case_15 output_15
	#perf stat ./exec_test test_case_20 output_20


gen_test: test_case_gen
	for i in {11..21}; do \
  		echo $$i; \
		./test_case_gen $$i ; \
		#               ./exec test_case_$$i out_$$i; \
	done;


#	./test_case_gen 10
#	./test_case_gen 15
#	./test_case_gen 20

test_case_gen: test_case_gen.o
	$(CC) $(CFLAGS) -o test_case_gen test_case_gen.o $(LDFLAGS)

test_case_gen.o:test_case_gen.cpp
	$(CC) $(CFLAGS) -c test_case_gen.cpp