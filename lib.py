### LIBRARY OF MODELS

from sequence import *
import numpy as np

N = -1 # placeholder value

mu = 1.0
JC69 = Model('JC69') # Jukes & Cantor 1969
JC69.addRM(np.array([[N, 1, 1, 1],
                     [1, N, 1, 1],
                     [1, 1, N, 1],
                     [1, 1, 1, N]])*(mu/4))

k = 2.0 # Transition-transversion ratio
K80 = Model('K80') # Kimura 1980
K80.addRM(np.array([[N, k, 1, 1],
                    [k, N, 1, 1],
                    [1, 1, N, k],
                    [1, 1, k, N]]))

models = [JC69, K80]

def simShell(seq: Sequence, tmax: float):
    seq1 = seq
    seq2 = seq.copy()
    all_data = {seq1: [seq1], seq2: [seq2]}
    for s in all_data.keys():
        t = 0
        s_sim = s
        while t < tmax:
            dt = np.random.exponential(-1/(seq.seq_l*seq.model.rate_mat[0,0]))
            t += dt
            s_sim = s_sim.sim(dt)
            all_data[s].append(s_sim)
    return all_data

def allGenDists(data: dict[Sequence, list[Sequence]]):
    cloned_data = {k: data[k].copy() for k in data.keys()}
    d_keys = list(data.keys())
    for i in [0,1]:
        other_s = data[d_keys[i-1]]
        for s in data[d_keys[i]][1:]:
            clone = Sequence()
            for o_s in other_s:
                if s.loc_t > o_s.loc_t:
                    clone = o_s.copy()
                    clone.loc_t = s.loc_t
            if clone.model is not None: cloned_data[d_keys[i-1]].append(clone)
        cloned_data[d_keys[i-1]].sort(key=lambda x: x.loc_t)
    seq1s, seq2s = cloned_data.values()
    ts = [s.loc_t for s in seq1s]
    ds = [seq1s[i].genDist(seq2s[i]) for i in range(len(seq1s))]
    return ts, ds