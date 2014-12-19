#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/Potts/TestPottsMutations.hpp
#

start_run=0
end_run=60

for (( i=start_run ; i<end_run; i = i+10))
do
    echo "Beginning run $i"
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
    # ">" directs std::cout to the file.
    # "2>&1" directs std::cerr to the same place.
    # "&" on the end lets the script carry on and not wait until this has finished.
    nice -20 ../../build/optimised/Potts/TestPottsCryptMutantRunner -run_index $i -num_runs 10 > output/run_${i}_output.txt 2>&1 &
done

echo "Jobs submitted"


