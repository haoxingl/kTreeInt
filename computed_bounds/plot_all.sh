#!/bin/bash

# Default arguments
bound_type="ours"
show=0

# Loop over values of m
for m in 64 128 256 512
do
  # Loop over values of prob
  for prob in 0.001 0.01 0.5 0.99
  do
    for zm in 0 1
    do
      for type in "type2_ub" "type2_lb"
      do
        python3 plot_bounds.py -t $type -zm $zm -m $m -k 4 8 16 32 64 128 256 512 1024 2048 4096 -prob $prob -b $bound_type -show $show
      done
    done
  done
done

# Loop over values of m
for m in 64 128 256 512
do
  # Loop over values of k
  for k in 4 8 16 32 64 128 256 512 1024
  do
    for zm in 0 1
    do
      python3 plot_bounds.py -t type1 -zm $zm -m $m -k $k -b $bound_type -show $show
    done
  done
done