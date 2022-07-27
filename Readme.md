# Readme

## Source code info
Programming Language: `C++` 

Additional Programming Language Info: 
- Compiler Info: `gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04)`  
- Packages/Libraries Needed: `CMake 2.8`, `makefile`

## Hardware Info
- Processor: 2  `Intel(R) Xeon(R) Silver 4110 CPU @ 2.10GHz` processors
- Caches: 3 level caches (`512KiB L1 cache`, `8MiB L2 cache`, `11MiB L3 cache`) for each processor
- Memory: `256GiB System Memory` (8 \* `32GiB DIMM DDR4 Synchronous 2666 MHz (0.4 ns)`)
- Secondary Storage: HDD, `6001GB TOSHIBA MG04ACA6`, (interface speed: `6.0 Gbit/s Max.`, rotation speed: `7,200 rpm `, average latency time: `4.17 ms `, buffer size: `128 MiB`, data transfer speed: `205 MiB/s `) write speed: `210-280 MiB/s`, read speed: `320-350 MiB/s`
- Network: there is no network usage in our experiments

## Dataset info
Most datasets used in our experiments are from the 
[KONECT](http://konect.cc/networks/),
[Network Repository](https://networkrepository.com),
and [SNAP](https://snap.stanford.edu/data/).

## How to run the code
After download and compile the code, the users can type `./LDScvx ` to see the usage of the program.
For example, `./LDScvx -g AM/graph.txt -k 5`  is running `LDScvx` on AM dataset to find top-5 LDS's.

## How to reproduce the experiments
To reproduce Figure 8, run the `LDScvx` (ours) and `LDSflow` (the codes for KDD15) for the 9 datasets. For example, `./LDScvx -g AM/graph.txt -k 5`  is running `LDScvx` on AM dataset to find top-5 LDS's. The last line of the output contains the running time.

For Figure 9, run the two algorithms for `k = 5, 10, 15, 20, 25`, respectively, and record the running time.

For Figure 10, we run `LDScvx` on 5 synthetic datasets and record the running time.

For Figure 11, we use `/usr/bin/time -v ./LDScvx -g AM/graph.txt -k 5` to catch the memory usage of `LDScvx` on AM dataset. Similarly, the memory usages of both algorithms over all datasets can be recorded. 

For Figure 12, run the `LDScvx` (ours) and `LDSflow` (the codes for KDD15) for the 9 datasets. The last 3 line contains the running time numbers of each part, which can be used to compute the proportion plot in Figure 12.

For Figure 13, run the `LDScvx`, `Greedy`, and `FDS` following the settings in the figure, and record the subgraph sizes, which are reported in the files named by each dataset (e.g., AM) and the value of `k` under the `output/` subdirectory.