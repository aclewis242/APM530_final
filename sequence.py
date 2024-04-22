### NUCLEOTIDE BASE SEQUENCE

from model import *
import numpy as np
import scipy as sp

A, C, G, T = 0, 1, 2, 3 # Nucleotide base indices
ACGT = [A, C, G, T]     # List of bases
ACGT_str = 'ACGT'       # String representation

class Sequence:
    '''
    The class representing a sequence of nucleotide bases (A, C, G, T).
    '''
    seq_l = -1                  # Length of base sequence
    model = None                # Which model to use
    base_freqs = np.empty(4)    # Initial base frequencies
    rate_mat = np.empty(1)      # Instantaneous rate matrix
    seq = []                    # Base sequence (as a list)
    loc_t = 0                   # Tracks total simulation time (local)

    def __init__(self, sl: int=None, mdl: Model=None, sq: list=[], lt: float=0) -> None:
        '''
        Builds the Sequence object.

        ### Parameters
        sl: Sequence length
        mdl: Substitution model
        sq: List of bases
        lt: Total local simulation time
        '''
        self.loc_t = lt
        if sl is not None: self.seq_l = sl
        if mdl is not None:
            self.model = mdl
            self.base_freqs = mdl.base_freqs
            self.rate_mat = mdl.rate_mat
            if len(sq): self.seq = sq
            else: self.gen(mdl.base_freqs)
    
    def gen(self, bfs: list) -> None:
        '''
        Generates a base sequence from the given set of frequencies.

        ### Parameters
        bfs: List of base frequencies
        '''
        self.seq = [np.random.choice(ACGT, p=bfs) for i in range(self.seq_l)]

    # TODO: matrix exponential is (probably) inefficient - find a way around using it if necessary
    def trans_ps(self, t: float) -> np.ndarray:
        '''
        Finds transition probability matrix from the rate matrix of the model at a given time.

        ### Parameters
        t: The length of time passed since the start
        '''
        return sp.linalg.expm(self.rate_mat*t)
    
    def sim(self, t: float) -> 'Sequence':
        '''
        Simulates the evolution of the sequence over a given time period.

        ### Parameters
        t: The length of the time period
        '''
        tpm = self.trans_ps(t)
        seq = [np.random.choice(ACGT, p=tpm[b]) for b in self.seq]
        return Sequence(sl=self.seq_l, mdl=self.model, sq=seq, lt=self.loc_t+t)
    
    def genDist(self, seq2: 'Sequence') -> float:
        '''
        Calculates observed distance (p-distance) between the current and original sequences. Serves as a measure of genetic drift.
        '''
        if self.seq_l != seq2.seq_l: print('Sequences are of different lengths!')
        return sum([self.seq[i] != seq2.seq[i] for i in range(self.seq_l)])/self.seq_l
    
    def copy(self) -> 'Sequence':
        return Sequence(sl=self.seq_l, mdl=self.model, sq=self.seq, lt=self.loc_t)
    
    @property
    def totals(self) -> list:
        '''
        A 4-element list describing the current macrostate of the sequence. (i.e., the total #s of each nucleotide)
        '''
        return [self.seq.count(b) for b in ACGT]

    def printTotals(self) -> None:
        '''
        Outputs the current macrostate (nucleotide #s) to the console.
        '''
        [print(f'{ACGT_str[b]}: {self.totals[b]}') for b in ACGT]

    def __str__(self) -> str:
        return ''.join([ACGT_str[b] for b in self.seq])
    
    def __repr__(self) -> str:
        return self.__str__()