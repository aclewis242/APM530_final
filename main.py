### MAIN FILE

from lib import *
import matplotlib.pyplot as plt

def runGDs(mdl: Model, seq_l: int=50, t_max: float=5, to_show: bool=True):
    '''
    Run the genetic drift simulation for a single model and plot the results.

    ### Parameters
    mdl: The model to use (refer to lib.py)
    seq_l: The length of the nucleotide base sequence
    t_max: The maximum amount of time to use in the simulation
    to_show: Whether or not to display the graph produced

    ### Returns
    ts: The times visited by the simulation
    ds: The genetic drift values at each time
    '''
    ## Simulation
    seq = Sequence(seq_l, mdl)
    ts, ds = allGenDists(simShell(seq, t_max))

    ## Plot results
    plt.plot(ts, ds)
    plt.title(rf'{mdl.name}: genetic drift over time')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$d$')
    plt.savefig(f'gendrift_{mdl.name}.png')
    if to_show: plt.show()
    plt.close()

    return ts, ds

def runGDAll(seq_l: int=50, t_max: float=5, to_show: bool=True):
    '''
    Run the genetic drift simulation for all models and plot the results together.
    Returns a list containing pairs of times and genetic drifts for each model.

    ### Parameters
    seq_l: The length of the nucleotide base sequence
    t_max: Maximum simulation time
    to_show: Whether or not to display the results
    '''
    all_dat = {m: runGDs(m, seq_l, t_max, to_show=False) for m in models}

    [plt.plot(all_dat[m][0], all_dat[m][1], label=m) for m in models]
    plt.title(r'All models: genetic drift over time')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$d$')
    plt.legend()
    plt.savefig('all_gendrift.png')
    if to_show: plt.show()
    plt.close()

    return list(all_dat.values())

def runGDAvgs(num_samp: int=10):
    '''
    Runs the whole simulation for the specified number of times and plots statistical data about the results.

    ### Parameters
    num_samp: The number of times to run the simulation
    '''
    all_ds = []
    for i in range(num_samp):
        all_ds.append([np.mean(ad[1]) for ad in runGDAll(to_show=False)])
    all_means = []
    for i in range(len(all_ds[0])):
        ds = [ad[i] for ad in all_ds]
        all_means.append([np.mean(ds), np.std(ds)])
    
    [plt.bar(str(models[i]), all_means[i][0], yerr=all_means[i][1]) for i in range(len(all_means))]
    plt.title(r'Average genetic drifts')
    plt.xlabel(r'Model')
    plt.ylabel(r'$d_{avg}$')
    plt.savefig(f'avggd_{num_samp}s.png')
    plt.show()

def runBCs(mdl: Model, seq_l: int=50, t_max: float=5, to_show: bool=True):
    '''
    Runs the simulation for base counts and plots them all over time.

    ### Parameters
    mdl: The model to use
    seq_l: The length of the sequence
    t_max: Maximum simulation time
    to_show: Whether or not to display the graph
    '''
    seq = Sequence(seq_l, mdl)
    all_seqs = simShell(seq, t_max)[seq]
    all_tots = [s.totals for s in all_seqs]
    all_ts = [s.loc_t for s in all_seqs]

    [plt.plot(all_ts, np.array(all_tots)[:,b], label=ACGT_str[b]) for b in ACGT]
    plt.title(rf'{mdl}: nucleotide base counts')
    plt.xlabel(r'Time')
    plt.ylabel(r'Base count')
    plt.legend()
    plt.savefig(f'nbcs_{mdl}.png')
    if to_show: plt.show()

if __name__ == '__main__':
    ### Un/comment whichever line is of interest. The ones that run simulations are marked with *s.

    ### Individual model simulations
    mth = runBCs # runGDs: genetic drift // runBCs: base counts
    [mth(m) for m in models] # *

    ### Compiled simulation results
    # runGDAll() # *

    ### Simulation statistics
    # runGDAvgs(11) # *