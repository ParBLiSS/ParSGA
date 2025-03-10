ParSGA
========================================================================

ParSGA (**Par**allel **S**equence to **G**raph **A**ligner) is the first work-optimal parallel algorithm designed for semiglobal sequence alignment to genome graphs, such as variation graphs and splicing graphs. ParSGA leverages OpenMP to achieve high-performance alignment on multi-core CPUs, ensuring efficient utilization of available processing power. For in-depth information on the algorithm and its performance, please refer to our paper [below](#publication).

## Dependencies

- [cmake](https://cmake.org) version >= 3.1
- A C++ compiler with c++14 support
- [Google Protobuf](https://github.com/protocolbuffers/protobuf) library, also available using [conda](https://anaconda.org/anaconda/protobuf)

## Download and compile

The repository and external submodules can be downloaded using the recursive clone.

```sh
git clone --recursive <GITHUB_URL>
```

Next, compile the code using cmake utility:

```sh
mkdir build_directory && cd build_directory
cmake <OPTIONS> ../ParSGA
make
```

OPTIONS: 
1. `-DPROTOBUF_DIR=<path>` should provide *absolute* path to installation directory of google protobuf library.
2. Cmake will automatically look for default C/C++ compilers. To modify the default selection if needed, users can set the two variables `-DCMAKE_CXX_COMPILER=<path to C++ compiler>` and `-DCMAKE_C_COMPILER=<path to C compiler>`. 

After the compilation completes, expect the executables `preprocessing`, `parallel_exe`, and `serial_exe` in your build\_directory. 

## Usage

* Preprocess a .vg format graph to construct Compressed Sparse Row (CSR) format graph and character graphs and component labels:
```sh
preprocessing graph.vg -q output_path_prefix
```

* Run the serial algorithm to align a set of query sequences against the graph:
```sh
serial_exe graph.vg -q input_path_prefix output_path_prefix read_name
```

* Run the parallel ParSGA algorithm to align a set of query sequences against the graph:
```sh
parallel_exe graph.vg -q input_path_prefix output_path_prefix read_name
```

**Preprocessing file format:** The CSRs of the graph and character graphs and their component labels are stored with suffixes appended to the output_path_prefix given to the preprocessing executable. This format is assumed in the input to the serial_exe and parallel_exe executables with the given input_path_prefix.

**Read file format:** The read file is assumed to be in the format of \<input_path_prefix\>-read\<read_name\>.csr


## <a name=“publication”></a>Publication

- **Aranya Banerjee, Daniel Gibney, Helen Xu, and Srinivas Aluru**. "[A Work-Optimal Parallel Algorithm for Aligning Sequences to Genome Graphs]()". *IEEE International Parallel and Distributed Processing Symposium (IPDPS) 2025*.