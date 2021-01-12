# A tool for BNSL benchmark test

Please use these codes for reproducing experiments.

## Installation

Python 
```bash
pip install networkx
pip install dwave-ocean-sdk
pip install pyarmor
```

## Usage

operation.py : solver of BNSL
```python
data = "asia" # dataset (asia, child, insurance, water, alarm, barley, hailfinder, hepar2) 
dir_ = "./" # directory
data = "asia" # dataset
iter_n = 10 # number of samplings
delta_ratio = 1e5 # penalty coefficient (delta = delta_lower * delta_ratio)
graphs = 5 # numbers of subgraphs for pseudo dataset (n = 5 * graphs)
Solver = "SA" # Annealing Solver (Simulated Annealing (SA), Quantum Annealer (QA), Digital Annealer (DA))
num_reads = 1000 # parameter for SA
annealing_time = 100 # parameter for QA
```