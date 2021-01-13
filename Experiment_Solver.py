dir_ = "./" # directory
data = "alarm" # dataset
N = 10000 # sample size
greedy = "true" # algorithm1
iter_n = 10 # number of samplings
Solver = "SA" # Annealing Solver (Simulated Annealing (SA), Quantum Annealer (QA))
delta_ratio = 1e5 # penalty coefficient (delta = delta_lower * delta_ratio)
num_reads = 1000 # parameter for SA
annealing_time = 100 # parameter for QA

import pickle
import pandas as pd
from statistics import mean, stdev
from EncodeSolver import Encode_Solver
score_list = []
with open(dir_ + data + "/" + data + "_" + str(N) + "_" + greedy + ".pickle", mode='rb') as f:
    H_list = pickle.load(f)
for num in range(iter_n):
    ES = Encode_Solver()
    ES.Hamiltonian(delta_ratio=delta_ratio, H=H_list[num])
    score_list += [ES.solver(Solver=Solver, num_reads=num_reads, annealing_time=annealing_time)]
if iter_n > 1:
    print("Result:")
    print("score:", "mean", mean(score_list), "std", stdev(score_list))
score_list = pd.DataFrame(score_list)
score_list.columns = ["score"]
score_list.to_csv(dir_ + data + "/" + data + "_10000_true_" + Solver + ".csv", index=False)