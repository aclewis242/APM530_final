### LIBRARY OF MODELS

from sequence import *
import numpy as np

N = -1 # placeholder value

m = 1.0
JC69 = Model('JC69') # Jukes & Cantor 1969
JC69.addRM(np.array([[N, 1, 1, 1],
                     [1, N, 1, 1],
                     [1, 1, N, 1],
                     [1, 1, 1, N]])*m)

k = 2.0 # Transition-transversion ratio
K80 = Model('K80') # Kimura 1980
K80.addRM(np.array([[N, 1, k, 1],
                    [1, N, 1, k],
                    [k, 1, N, 1],
                    [1, k, 1, N]])*m)

bfs = {A: 0.1, C: 0.5, G: 0.2} # Base frequencies (T implicit, as they must sum to 1)
F81 = Model('F81', bfs) # Felsenstein 1981
F81.addRM(JC69.rate_mat*4)

HKY85 = Model('HKY85', bfs) # Hasegawa, Kishino, Yano 1985
HKY85.addRM(K80.rate_mat*4)

g = 0.5 # pyrimidine/purine transition ratio
t = 2*k/(1+g) # net purine coefficient (g*t = pyrimidine)
TN93 = Model('TN93', bfs) # Tamura & Nei 1993
TN93.addRM(np.array([[N, 1, t, 1],
                     [1, N, 1, t*g],
                     [t, 1, N, 1],
                     [1, t*g, 1, N]])*m)

[a, b, c, d, e, f] = [0.7, 1.3, 1.1, 0.9, 1.5, 0.6] # coefficients for GTR
GTR = Model('GTR', bfs) # Tavaré 1986
GTR.addRM(np.array([[N, a, b, c],
                    [a, N, d, e],
                    [b, d, N, f],
                    [c, e, f, N]]))

models = [JC69, K80, F81, HKY85, TN93, GTR]

def simShell(seq: Sequence, tmax: float):
    '''
    The core of all simulation methods. Models the evolution of two sequences from the given ancestral sequence over time.

    ### Parameters
    seq: The initial ancestral sequence
    tmax: Maximum simulation time
    '''
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
    '''
    Takes the two sequences produced by the simulation shell and turns them into a set of genetic distances.

    ### Parameters
    data: A dict linking the two sequence objects to the lists of them over time.
    (Note: This should almost always be the output of simShell()!)

    ### Returns
    ts: The list of times visited by both simulations
    ds: The genetic distance at each time
    '''
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
    ps = [seq1s[i].genDist(seq2s[i]) for i in range(len(seq1s))]
    ds = [ps[0]]
    [ds.append(ps[i]+ds[i-1]) for i in range(1,len(ps))]
    return ts, ds