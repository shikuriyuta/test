data = "hailfinder" # dataset
iter_n = 10 # number of samplings

import pickle
with open(data + "_10000_true.pickle", mode='rb') as f:
    H_list = pickle.load(f)
for num in range(iter_n):
    X = H_list[num][0] # variable names
    U_list = H_list[num][1] # candidate parent sets
    V_list = H_list[num][2] # candidate parent sets
    s0_list = H_list[num][3] # cofficient for score component of Hamiltonian
    s1_list = H_list[num][4] # cofficient for score component of Hamiltonian
    s2_list = H_list[num][5] # cofficient for score component of Hamiltonian
    s12_list = H_list[num][6] # cofficient for score component of Hamiltonian
    delta_lower = H_list[num][7] # lower boundary of delta