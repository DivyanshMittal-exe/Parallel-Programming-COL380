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






