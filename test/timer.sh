#!/bin/bash

# Initialize a variable to store the total execution time
total=0

# Loop 100 times
for i in {1..100}
do
    # Get the current time in seconds (with nanoseconds for precision) and store it in 'start'
    start=$(date +%s.%N)

    # Run the program with input redirected from 'n40Ca.in'
    ./COLOSS < d93Nb.in 

    # Get the current time in seconds again and store it in 'end'
    end=$(date +%s.%N)

    # Calculate the difference between 'end' and 'start' (this is the execution time of the program)
    # and store it in 'diff'
    diff=$(echo "$end - $start" | bc)

    # Add 'diff' to 'total'
    total=$(echo "$total + $diff" | bc)
done

# Print the total time taken
echo "Total time taken: $total seconds"
