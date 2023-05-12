Explanation of files:
Makefiles has compilation for each file. Edit them when compiling based on compiler flag you decided to use.

# Syntethic Data Processing : 
created for debug purposes. each generates vectors as input and works on them. all variations below uses openmp rather than stapl.
- find_bug : the version with static scheduling and chunk size. implements load balance of 2 variations, either in parallel or as preprocessing. commits must have both versions, current one is parallel lb.
- find_bug_schedule:  the version with dynamic scheduling and chunk size
- find_debug-serial : completely serial version without openmp.

# Real Data Processing : 
uses real data, either mock curves (transformed to fits in google colab for using the inherent data cleaning functions of fits library) or real fits files downloaded from kepler website. data reader changes based on the type. code has comments for it.
below files should be almost identical to above files in the sense that changes must be applied to them after debug in above files.
-  ompBLS.cpp : implemented with openmp 
-  omp-serialBLS : serial version of the above code
-  parallelBLS : uses stapl. run with mpi as the code uses mpi wall time to measure execution time
-  serial.cpp :  oldest serial version of the above with no optimizations. direct c++ version for naive python. use above serial codes for speedup comparisons.

# Experiments: 
Logs folder is old and used for solving a bug, ignore it. \
Experiments are under experiments directory. Their more sophisticated categorization is done in Campus Cluster rather than Parasol. \
Campus cluster experiments (where we achieve good  speedup with a max of 14 for 16 threads) is under cc-experiments directory. \
CC organizes experiments based on the input size, in each folder there are different runs for variety of implementations (scheduling, load balancing)

# Notes:
STAPL Solutions: \
\
Working STAPL solution can be found in commit with id : commit de2e9c2. See parallelBls.cpp and it has balanced distribution. \
Stapl balanced distribution is as slow as no distirbution specified. \
The most plain  and working STAPL solution can be found in commit with id : commit 0d01ca3. Has no distribution specified Vectors are initialized locally but automatically distributed by STAPL runtime system. Slow due to contention.

Performance Observations: 


for size=873, with both syntethic and mock file, debug and ompBLS programs finish around 6*10^-3 seconds with 4 threads.
STAPL version parallelBLS with 4 MPI processes take around 3.6 seconds for the same mock file of 873 size when using stapl arrays.  for this one bottleneck is not reduction but the communications, as stapl vectors are not maintained local to process. removing reduction over custom defined MPI reduction only faster by 0.002 second.
Stapl local vectors by using vector add is as slow due to incurred communication.
With std::vectors however, 4 MPI processes finish 0.0006 seconds for the same mock file. The same version gives stapl runtime errors with big datafiles such as 112k.
All these versions are committed to parallelBLS.cpp file separately.

Implementing the load balance is no trivial as map functions requires wrapping local data of processes into nested vectors before distributing them. Future library implementation is the last work that is started for the sake of overlapping load communication and computation. \
Based on trials, map function unrolls arguments of stapl vectors of std vectors, but communciations via future cannot be performed due to std. \
Alternatively, passing stapl array views of stapl views of data does not unroll to stapl views of data in each location.