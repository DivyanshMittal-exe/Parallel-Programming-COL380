#!/bin/bash

# Create a file to store the results
results_file="results.txt"
echo "Size Cores Time_Elapsed Non_zero_blocks" > $results_file

# Loop through the sizes
#!/bin/bash

# loop through size from 10 to 18
for i in {10..25}; do
  # Generate test case with size i
  ./test_case_gen $i

  # loop through number of cores from 1, 2, 4, 8, 16
  for j in {2,4,8,16}; do
    # set the number of threads for OpenMP
#    export OMP_NUM_THREADS=$j

    TIMEFORMAT=%R

#    time_elapsed=$(time OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm) 2>&1

#    time_elapsed=$( TIMEFORMAT="%R"; { time OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm; } 2>&1 )
#    time_elapsed=$( { OMP_NUM_THREADS=$j perf stat ./exec test_case_$i output_mmmm; } 2>&1 | grep user | awk '{print $1}' )

    time_elapsed=$( { OMP_NUM_THREADS=$j perf stat ./exec test_case_$i output_mmmm; } 2>&1 | grep -E 'user|sys' | awk '{print $1}' )


    # get time elapsed using perf stat
#    time_elapsed=$(perf stat -r 1 -o perf.txt OMP_NUM_THREADS=$j ./exec test_case_$i output_mmmm 2>&1 | awk '/time elapsed/ {print $1}' | awk '{print $1}')

    # run the comp file and get the output
    non_zero_blocks=$(./comp output_mmmm)

    # write the results to the file
    echo $i $j $time_elapsed $non_zero_blocks >> results.txt
    echo $i $j $time_elapsed $non_zero_blocks

  done

  # remove the generated test case
  rm test_case_$i
done

