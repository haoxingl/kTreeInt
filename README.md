# On Wagner's k-Tree Algorithm Over Integers

This repository contains the implementation of the computed bounds of success probability and the complexity of Wagner's k-Tree algorithm over integers. The results are presented in the paper "On Wagner's k-Tree Algorithm Over Integers" by Haoxing LIN and Prashant Vasudevan: https://arxiv.org/abs/2410.06856.

We also provide the implementation of evaluation code for the k-Tree algorithm, which compares its actual performance with our computed bounds.

## Prerequisites

- Python 3 (for computed bounds and plotting results): numpy, pandas, matplotlib
- CMake (version 3.10 or higher)
- C++ compiler with C++17 support (e.g., GCC or Clang)

## Parameters of the k-Tree Algorithm

The k-Tree algorithm has the following parameters (refer to the paper for more details):

- `m`: The range for the integers (m > 1)
- `k`: The number of input lists in the k-Tree algorithm (k > 1)
- `n`: The number of elements in each input list (n > 1)

## Project Structure

The computed bounds are implemented in the 'computed_bounds' directory, with the underlying functions in 'computed_bounds/proof_functions.py' and the plotting script in 'computed_bounds/plot_bounds.py'. There is also a bash script 'plot_all.sh' that allows you to plot multiple computed bounds with your specified parameters. The algorithms will output the computed bounds as pdf figures in the same directory.

Below are some examples of running the plotting script:

   ```
   python3 plot_bounds.py -t type1 -zm 0 -m 256 -k 512
   
   python3 python3 plot_bounds.py -t type2_ub -zm 0 -m 256 -k 4 8 16 32 64 128 256 512 1024 2048 4096 -prob 0.5

   python3 python3 plot_bounds.py -t type2_ub -zm 0 -m 256 -k 4 8 16 32 64 128 256 512 1024 2048 4096 -prob 0.5
   ```

The arguments are as follows:

- `-t`: The type of the bound. Choose from `type1`, `type2_ub`, and `type2_lb`. 'type1' is the bound for the success probability of the k-Tree algorithm, while 'type2_ub' and 'type2_lb' are the bounds for the sufficient and necessary complexity of the k-Tree algorithm to succeed.
- `-zm`: The boolean value of whether the bound is over $\mathbb{Z}_m$. Choose from 0 or 1.
- `-m`: The value of the exponent of `m`.
- `-k`: The value of `k`. For `type2_ub` and `type2_lb`, you can specify multiple values of `k` to plot the bounds for different values of `k`.

<!-- ## Building the Project

To build the project, follow these steps:

1. Open a terminal and navigate to the project root directory.

2. Create and build the project using CMake:

   ```
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   ```

   This will create the main executable in the `build` directory.

## Build Types

- Use `-DCMAKE_BUILD_TYPE=Debug` for debugging (adds debug symbols)
- Use `-DCMAKE_BUILD_TYPE=Release` for optimized build (applies -O3 optimization)

## Additional Notes

- The project is set up to use C++17 standard.
- For any changes to the build configuration, modify the `CMakeLists.txt` file in the project root. -->