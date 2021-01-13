# Efficient Conversion of Bayesian Network Structure Learning into Quadratic Unconstrainded Binary Optimization

Please use these codes for reproducing experiments.

## Installation

Python 
```bash
pip install bnlearn
pip install networkx
pip install pgmpy
pip install dwave-ocean-sdk
```

R 
```bash
install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")
```

Julia 
```bash
add Statistics 
add DataFrames 
add Feather 
add CSV
add BenchmarkTools
add Distributed
add Combinatorics
add SpecialFunctions
```

## Usage

Experiment_Dataset.r : sampling of benchmark data (asia, child, insurance, water, alarm, barley, hailfinder, hepar2)
```r
dir_ = "./" # directory
data = "asia" # dataset 
N = 10000 # sample size
iter_n = 10 # number of samplings
```

Experiment_Conversion.jl : conversion of BNSL (asia, child, insurance, water, alarm, barley, hailfinder, hepar2)
```julia
addprocs(10) # number of process
@everywhere dir_ = "./" # directory
data = "asia" # dataset
m = 2 # maximum size of parent set
N = 100 # sample size
greedy = false # algorithm1
iter_n = 10 # number of samplings
```

Experiment_Encode.py : encode of BNSL
```python   
dir_ = "./" # directory
data = "alarm" # dataset
N = 10000 # sample size
greedy = "true" # algorithm1
iter_n = 10 # number of samplings
```

Experiment_Solver.py : solver of BNSL
```python
dir_ = "./" # directory
data = "asia" # dataset
N = 10000 # sample size
greedy = "true" # algorithm1
iter_n = 10 # number of samplings
Solver = "SA" # Annealing Solver (Simulated Annealing (SA), Quantum Annealer (QA))
delta_ratio = 1e5 # penalty coefficient (delta = delta_lower * delta_ratio)
num_reads = 1000 # parameter for SA
annealing_time = 100 # parameter for QA
```