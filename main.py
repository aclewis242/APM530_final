### MAIN FILE

from sequence import *
from lib import *
import matplotlib.pyplot as plt

def run(mdl: Model=JC69, seq_l: int=1000, t_max: int=10):
    '''
    Run the simulation and plot the results.

    ### Parameters
    mdl: The model to use (refer to lib.py).
    seq_l: The length of the nucleotide base sequence.
    t_max: The maximum amount of time to use in the simulation.
    '''
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
    plt.savefig(f'seqev_{mdl}.png')
    plt.show()

if __name__ == '__main__':
    run(K80)