for i in {11..21}; do
    echo $i;
  ./test_case_gen $i ;
  ./exec test_case_$i out_$i

done;