Number of approaches:

1:Idea: Basic idea to make complete matrix and normal multiply
1:Why attempted: To check proper input output/processing for easy debug later
1:Results (Specify speedup): Too slow for larger input2, just used to make sure data processing works
1:Drawback: Too slow and does not use sparseness

2:Idea: Create a struct chunk with int x,y and data in a vector
2:Why attempted: To make use of sparseness
2:Results (Specify speedup): Much faster maked it easy to process the input2
2:Drawback: Was copying data and still io was bottleneck

3:Idea: Try use parallel execution policy for sort and async io
3:Why attempted: Speed up output and any sequential part left
3:Results (Specify speedup): Was unable to implement on CSS
3:Drawback: Execution policy didn't work on CSS and async io also gave problems so gave up

4:Idea: Now I switched to array as size fixed. When getting x,y from y,x I don't copy data now, also use template
4:Why attempted: Not copy data and speed up output
4:Results (Specify speedup): 1.7 times speed up, now compute final output in parallel, and just print sequentially fast
4:Drawback: Still major time is taken by output


Final Approach:
  Read chunk data in chunks and do memory dump which is fast. Create the transpose pair of chunks without copying data. Sort indices for easy lookup and calculate all possible pairs in the final output. Each row is a task and thus we dont have a race condition. Using sorted data allows faster indices lookup using binary search, thus it is easy to calculate all the possible (i,j) pairs for the 'i'th row. These final outputs are now converted to unsigned shorts parallely so that output can again be a simple memory dump.


Scalability analysis:

Size,Non_Zero_Blocks,Cores: 1,Cores: 2,Cores: 4,Cores: 8,Cores: 16
10,2948,990.104,850.089,630.065,500.052,740.076
11,5814,1300.132,1210.124,1060.108,1230.127,1170.121
12,11726,2650.268,2350.238,2200.223,2050.207,2050.208
13,23560,4190.423,4190.424,4370.441,4470.451,4370.446
14,46862,9300.936,10611.072,9760.987,8760.885,9170.929
15,94099,20452.066,17871.808,17101.724,18761.890,16971.712
16,187252,36473.664,38123.833,35273.551,35163.539,38963.916
17,375569,73927.429,75657.605,71677.204,68466.883,72657.304
18,750642,1379713.865,1239312.467,1204712.124,1143011.504,1290712.982

Size,Non_Zero_Blocks,Cores: 1,Cores: 2,Cores: 4,Cores: 8,Cores: 16
11,5855,280.032,190.021,200.023,150.017,150.017
12,11859,390.041,330.036,290.032,280.031,290.032
13,23354,780.082,660.069,570.061,550.059,560.060
14,46693,1550.160,1300.136,1200.126,1170.123,1110.117
15,93586,8350.846,2650.276,2610.273,2210.232,2450.258
16,187833,6560.674,5120.533,4710.491,4460.468,4670.488
17,375375,18491.880,10781.113,9460.981,8950.932,9400.976
18,751457,32593.321,21792.248,19502.021,18581.929,18381.907
19,1502294,56165.741,45384.676,40344.177,38063.948,37633.901

We see that we get best performance at 8 cores, because of the simple reason that CSS has 8 threads that it allots for running the code. Above that we need to incur overhead of thread switching


