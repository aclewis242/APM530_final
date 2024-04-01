### SEQUENCE EVOLUTION MODEL

import numpy as np

class Model:
    '''
    A container for the various pieces of information relevant to a substitution model.
    '''
    rate_mat = np.empty(1)      # Instantaneous rate matrix
    base_freqs = np.ones(4)/4   # Default initial base frequency is 'all equal'
    rps = {}                    # Dictionary of rate parameters
    name = 'blank'

    def __init__(self, nm: str, **kwargs: float) -> None:
        '''
        Builds the model container.

        ### Parameters
        nm: The name of the model
        (kwargs: any relevant rate parameters - optional)
        '''
        self.name = nm
        self.rps = kwargs
    
    def addRM(self, rm: np.ndarray) -> None:
        '''
        Adapts the input to a rate matrix.

        ### Parameters
        rm: The rate matrix (important: the diagonal values should be negative!)
        '''
        self.rate_mat = rm
        for i in range(4): self.rate_mat[i,i] = -rm[i,np.where(rm[i] > 0)[0]].sum()
    
    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return self.__str__()