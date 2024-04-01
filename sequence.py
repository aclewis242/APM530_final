### NUCLEOTIDE BASE SEQUENCE

import numpy as np
import scipy as sp

A, C, G, T = 1, 2, 3, 4 # Nucleotide base indices
ACGT = [A, C, G, T]     # List of bases

class Sequence:
    seq_l = -1                  # Length of base sequence
    rate_mat = np.empty(4,4)
    base_rates = np.empty(4)
    seq = []

    def __init__(self, sl, rm, brs):
        self.seq_l = sl
        self.rate_mat = rm
        self.base_rates = brs
        self.seq = [np.random.choice(ACGT, p=brs) for i in range(sl)]