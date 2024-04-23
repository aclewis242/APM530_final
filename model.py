### SEQUENCE EVOLUTION MODEL

import numpy as np

class Model:
    '''
    A container for the various pieces of information relevant to a substitution model.
    '''
    rate_mat = np.array([])     # Instantaneous rate matrix
    base_freqs = np.ones(4)/4   # Default initial base frequency is 'all equal'
    name = 'blank'

    def __init__(self, nm: str, bfs: dict=None):
        '''
        Builds the model container.

        ### Parameters
        nm: The name of the model
        bfs: A dict linking base ID to the corresponding frequency
        '''
        self.name = nm
        if bfs is not None:
            bfs_temp = bfs.copy()
            bfs_temp[6-sum(bfs.keys())] = 1 - sum(bfs.values())
            self.base_freqs = np.array([bfs_temp[i] for i in bfs_temp])
    
    def addRM(self, rm: np.ndarray):
        '''
        Adapts the input to a rate matrix.

        ### Parameters
        rm: The rate matrix (important: the diagonal values should be negative!)
        '''
        self.rate_mat = rm
        for i in range(4):
            self.rate_mat[i,np.where(rm[i] > 0)[0]] *= self.base_freqs[np.where(rm[i] > 0)[0]]
            self.rate_mat[i,i] = -rm[i,np.where(rm[i] > 0)[0]].sum()

    def __str__(self):
        '''
        String representation of the model (its name).
        '''
        return self.name
    
    def __repr__(self):
        '''
        Printable representation of the model (its name).
        '''
        return self.__str__()