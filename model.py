### SEQUENCE EVOLUTION MODEL

import numpy as np

class Model:
    '''
    A container for the various pieces of information relevant to a substitution model.
    '''
    rate_mat = np.empty(1)      # Instantaneous rate matrix
    base_freqs = np.ones(4)/4   # Default initial base frequency is 'all equal'
    tv_ratio = 1                # Transition/transversion rate ratio
    name = 'blank'

    def __init__(self, nm: str, tr=None) -> None:
        '''
        Builds the model container.

        ### Parameters
        nm: The name of the model
        tr: Transition/transversion rate ratio
        '''
        self.name = nm
        if tr is not None: self.tv_ratio = tr
    
    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return self.__str__()