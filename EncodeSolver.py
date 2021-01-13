import numpy as np
import pylab
import bnlearn as bn
import networkx as nx
from pgmpy.estimators import BDeuScore
import neal
import dimod
from dimod.reference.samplers import ExactSolver
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
from dwave_qbsolv import QBSolv
#from dwave.cloud import Client

class Encode_Solver:
    
    def __init__(self, D=None, G0=None):
        self.D = D
        self.G0 = G0
        self.score0 = 0.0
    
    def encode(self, EC, ESS=1.0):
        n = int((len(EC.columns) - 2) / 3)        
        self.X = EC.columns[(2 * n + 1):(3 * n + 1)]
        self.U_list = [] 
        self.V_list = [] 
        self.s0_list = [] 
        self.s1_list = [] 
        self.s2_list = [] 
        self.s12_list = [] 
        self.delta_lower = 0.0 
        self.score0 = 0.0
        count = 0
        if self.D is not None:
            BDeu = BDeuScore(self.D, equivalent_sample_size = ESS)
        for (i, X_) in enumerate(self.X):
            ES_ = EC[EC["i"] == i + 1]            
            UV = [set([self.X[j] for j in range(n) if ES_.iloc[h, 1 + 2 * n + j] == 1]) for h in range(len(ES_))]
            foo = []
            U = [U_ for U_ in [set([self.X[j] for j in range(n) if ES_.iloc[h, 1 + j] == 1]) for h in range(len(ES_))] 
                                                          if (U_ not in foo and not foo.append(U_)) and (U_ != set())]
            foo = []
            V = [V_ for V_ in [set([self.X[j] for j in range(n) if ES_.iloc[h, 1 + n + j] == 1]) for h in range(len(ES_))] 
                                                              if (V_ not in foo and not foo.append(V_)) and (V_ != set())]
            for (h, UV_) in enumerate(UV):
                if UV_ == set():
                    s0 = ES_.iloc[h, -1]
                    break
            s1 = []
            for U_ in U:
                for (h, UV_) in enumerate(UV):
                    if UV_ == U_:
                        s1 += [ES_.iloc[h, -1] - s0]
                        break
            s2 = []
            for V_ in V:
                for (h, UV_) in enumerate(UV):
                    if UV_ == V_:
                        s2 += [ES_.iloc[h, -1] - s0]
                        break
            s12 = []
            for (h1, U_) in enumerate(U):
                s12_ = []
                for (h2, V_) in enumerate(V):
                    for (h, UV_) in enumerate(UV):
                        if UV_ == (U_ | V_):
                            s12_ += [ES_.iloc[h, -1] - s1[h1] - s2[h2] - s0]
                            break
                s12 += [s12_]
            self.U_list += [U]
            self.V_list += [V]
            self.s0_list += [s0]
            self.s1_list += [s1]
            self.s2_list += [s2]
            self.s12_list += [s12]
            if s1 == [] and s2 == [] :
                self.delta_lower = self.delta_lower
            elif s1 == []:
                self.delta_lower = max(self.delta_lower, -min(s2))
            elif s2 == []:
                self.delta_lower = max(self.delta_lower, -min(s1))
            else:
                self.delta_lower = max(self.delta_lower, 
                                   max([-s1[i]-min(np.array(s12)[i,:]) for i in range(len(s1))]),
                                   max([-s2[j]-min(np.array(s12)[:,j]) for j in range(len(s2))]))
            if self.G0 is not None:
                ad = self.G0['adjmat'][X_]
                UV0 = set([ad.index[j] for j in range(n) if ad[j] == True])
                for UV_ in UV:
                    if UV_ == UV0:
                        count += 1
                        break
                if self.D is not None:
                    self.score0 += - BDeu.local_score(X_, list(UV0))
        if self.G0 is not None:
            print("percentage of true parent sets:", count / n)
            if self.D is not None:
                print('score0:', self.score0)
        return [self.X, self.U_list, self.V_list, self.s0_list, self.s1_list, self.s2_list, self.s12_list, self.delta_lower]
    
    def Hamiltonian(self, delta_ratio=2.0, H=None):
        if H is not None:
            self.X = H[0]
            self.U_list = H[1]
            self.V_list = H[2] 
            self.s0_list = H[3] 
            self.s1_list = H[4] 
            self.s2_list = H[5] 
            self.s12_list = H[6] 
            self.delta_lower = H[7]
        n = len(self.X)
        delta_lower_ = self.delta_lower * delta_ratio
        delta1 = delta_lower_ 
        delta2 = delta_lower_ * 3 
        delta3 = delta2 * (n - 2) / 3 
        self.s0_sum = sum(self.s0_list)
        self.H = dict()
        for i in range(n):
            if self.U_list[i] != []:
                for ii in range(len(self.U_list[i])):
                    if i > 0:
                        c = sum([1 for iii in range(i) if self.X[iii] in self.U_list[i][ii]])
                    else:
                        c = 0
                    self.H.update({((1,ii,i),(1,ii,i)): self.s1_list[i][ii] + delta3 * c})
            if self.V_list[i] != []:
                for ii in range(len(self.V_list[i])):
                    if i > 0:
                        c = sum([1 for iii in range(i) if self.X[iii] in self.V_list[i][ii]])
                    else:
                        c = 0
                    self.H.update({((2,ii,i),(2,ii,i)): self.s2_list[i][ii] + delta3 * c})
            if self.U_list[i] != [] and self.V_list[i] != []:
                for ii in range(len(self.U_list[i])):
                    for iii in range(len(self.V_list[i])):
                        self.H.update({((1,ii,i),(2,iii,i)): self.s12_list[i][ii][iii]})
            if len(self.U_list[i]) >= 2:
                for ii in range(len(self.U_list[i])-1):
                    for iii in range(ii+1,len(self.U_list[i])):
                        self.H.update({((1,ii,i), (1,iii,i)): delta1})
            if len(self.V_list[i]) >= 2:
                for ii in range(len(self.V_list[i])-1):
                    for iii in range(ii+1,len(self.V_list[i])):
                        self.H.update({((2,ii,i), (2,iii,i)): delta1})
        penalty = np.zeros((n, n))
        for i in range(n-2):
            for ii in range(i+1,n-1):
                for iii in range(ii+1,n):
                    penalty[i, iii] += delta2
                    self.H.update({((i,iii), (ii,iii)): -delta2})
                    self.H.update({((i,iii), (i,ii)): -delta2})
                    self.H.update({((i,ii), (ii,iii)): delta2})
        for i in range(n-2):
            for iii in range(i+2,n):
                self.H.update({((i,iii), (i,iii)): penalty[i, iii]})
        for i in range(1,n):
            for ii in range(i):
                if self.U_list[ii] != []:
                    for iii in range(len(self.U_list[ii])):
                        if self.X[i] in self.U_list[ii][iii]:
                            self.H.update({((1,iii,ii),(ii,i)): delta3})
                if self.U_list[i] != []:
                    for iii in range(len(self.U_list[i])):
                        if self.X[ii] in self.U_list[i][iii]:
                            self.H.update({((1,iii,i),(ii,i)): - delta3})
                if self.V_list[ii] != []:
                    for iii in range(len(self.V_list[ii])):
                        if self.X[i] in self.V_list[ii][iii]:
                            self.H.update({((2,iii,ii),(ii,i)): delta3})
                if self.V_list[i] != []:
                    for iii in range(len(self.V_list[i])):
                        if self.X[ii] in self.V_list[i][iii]:
                            self.H.update({((2,iii,i),(ii,i)): - delta3})
        return self.H
    
    def solver(self, Solver="SA", num_reads=1000, annealing_time=100, display=False, compare=False):
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.H, self.s0_sum)
        if Solver == "Exact":
            self.response = ExactSolver().sample(bqm)
        elif Solver == "SA":
            self.response = neal.SimulatedAnnealingSampler().sample(bqm, num_reads=num_reads)
        elif Solver == "QBSolv":
            self.response = QBSolv().sample(bqm)
        elif Solver == "QA":
            self.response = EmbeddingComposite(DWaveSampler()).sample(bqm, annealing_time=annealing_time)
        sample = None
        for a, r in enumerate(self.response):
            if a == 0:
                sample = r
                break
        edges = []
        fuga = 0
        for i, X_ in enumerate(self.X):
            if self.U_list[i] != []:
                hoge = 0
                for ii in range(len(self.U_list[i])):
                    if sample[(1,ii,i)] == 1:
                        for iii in list(self.U_list[i][ii]):
                            edges += [(iii, X_)]
                        hoge += 1
                if hoge > 1:
                    fuga += 1
            if self.V_list[i] != []:
                hoge = 0
                for ii in range(len(self.V_list[i])):
                    if sample[(2,ii,i)] == 1:
                        for iii in list(self.V_list[i][ii]):                            
                            edges += [(iii, X_)]
                        hoge += 1
                if hoge > 1:
                    fuga += 1
        G = nx.DiGraph(edges)
        for X_ in self.X:
            G.add_node(X_)
        if fuga > 0:
            print("Constraints Error: p is more than 2.")
        if nx.is_directed_acyclic_graph(G) == False:            
            print("Constraints Error: Not DAG")
        score = self.s0_sum
        for i, X_ in enumerate(self.X):
            UV = set()
            for e in edges:
                if e[1] == X_:
                    UV = UV | {e[0]}
            Z = set()
            for U in self.U_list[i]:
                Z = Z | U
            Z_ = set()
            for V in self.V_list[i]:
                Z_ = Z_ | V
            for ii, U in enumerate(self.U_list[i]):
                if U == UV & Z:
                    score += self.s1_list[i][ii]
            for iii, V in enumerate(self.V_list[i]):
                if V == UV & Z_:
                    score += self.s2_list[i][iii]
            for ii, U in enumerate(self.U_list[i]):
                for iii, V in enumerate(self.V_list[i]):
                    if U | V == UV:
                        score += self.s12_list[i][ii][iii]
        print('score:', score)
        if display != False:
            pylab.figure()
            pos = nx.spring_layout(G)
            nx.draw_networkx_nodes(G, pos, node_size=8, node_color="w")
            nx.draw_networkx_edges(G, pos, width=1)
            nx.draw_networkx_labels(G, pos, font_size=8, font_color="r")
            pylab.show()
        G = bn.make_DAG(list(G.edges))
        if (self.G0 is not None) and (compare == True):
            bn.compare_networks(G, self.G0, showfig=False)
        return score