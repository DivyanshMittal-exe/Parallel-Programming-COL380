#!/bin/bash

results_file="results3.txt"
echo "Size Cores Time_Elapsed Non_zero_blocks" > $results_file


for i in {11..25}; do
  ./test_case_gen $i

  for j in {1,2,4,8,16}; do

    TIMEFORMAT=%R

#    time_elapsed=$(time OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm) 2>&1

    time_elapsed=$( TIMEFORMAT="%R"; { time OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm; } 2>&1 )
  #  time_elapsed=$( { OMP_NUM_THREADS=$j perf stat ./exec test_case_$i output_mmmm; } 2>&1 | grep user | awk '{print $1}' )

    # time_elapsed=$( { OMP_NUM_THREADS=$j perf stat ./exec test_case_$i output_mmmm; } 2>&1 | grep -E 'user|sys' | awk '{print $1}' )
#    time_elapsed=$(OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm)


    non_zero_blocks=$(./comp output_mmmm)

    echo $i $j $time_elapsed $non_zero_blocks >> $results_file
    echo $i $j $time_elapsed $non_zero_blocks

  done

  rm test_case_$i
done

