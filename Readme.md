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