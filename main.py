### MAIN FILE

from sequence import *
from lib import *
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ## User-specified parameters
    seq_l = 1000    # Length of nucleotide base sequence
    mdl = JC69      # Model to use (refer to lib.py)
    t_max = 10      # Max time

    ## Simulation
    seq = Sequence(seq_l, mdl)
    seq.gen([1, 0, 0, 0]) # To better observe evolution towards equilibrium
    print('\tt = 0')
    seq.printTotals()
    all_seqs = [seq.sim(t) for t in range(t_max)]
    print(f'\tt = {t_max}')
    all_seqs[-1].printTotals()
    all_tots = np.array([s.totals for s in all_seqs])

    ## Plot results
    [plt.plot(range(t_max), all_tots[:,b], label=ACGT_str[b]) for b in ACGT]
    plt.title(f'Sequence evolution: {mdl}')
    plt.legend()
    plt.xlabel(r'$t$')
    plt.ylabel('Total base count')
    plt.show()