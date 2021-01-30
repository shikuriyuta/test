dir_ = "./" # directory
data = "water" # dataset
N = 10000 # sample size
greedy = "false" # algorithm1
iter_n = 10 # number of samplings

import pickle
import pandas as pd
from EncodeSolver import Encode_Solver
H_list = []
for num in range(iter_n):
    ES = Encode_Solver()
    EC = pd.read_csv(dir_ + data + "/" + data + "_" + str(N) + "_" + greedy + "_" + str(num + 1) + ".csv")
    H = ES.encode(EC)
    H_list += [H] 
with open(dir_ + data + "/" + data + "_" + str(N) + "_" + greedy + ".pickle", mode='wb') as f:
    pickle.dump(H_list, f)