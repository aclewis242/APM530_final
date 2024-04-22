### MAIN FILE

from lib import *
import matplotlib.pyplot as plt

def run(mdl: Model=JC69, seq_l: int=50, t_max: float=10):
    '''
    Run the simulation and plot the results.

    ### Parameters
    mdl: The model to use (refer to lib.py).
    seq_l: The length of the nucleotide base sequence.
    t_max: The maximum amount of time to use in the simulation.
    '''
    ## Simulation
    seq = Sequence(seq_l, mdl)
    ts, ds = allGenDists(simShell(seq, t_max))

    ## Plot results
    plt.plot(ts, ds)
    plt.title(f'{mdl.name}: genetic drift over time')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$d$')
    plt.savefig(f'gendrift_{mdl.name}.png')
    plt.show()

if __name__ == '__main__':
    [run(m) for m in models]
    # seq = Sequence(10, JC69)
    # seq.gen([1, 0, 0, 0]) # To better observe evolution towards equilibrium
    # print(seq.sim(1))