### LIBRARY OF MODELS

from model import *
import numpy as np

N = -1 # placeholder value

mu = 1
JC69 = Model('JC69') # Jukes & Cantor 1969
JC69.addRM(np.array([[N, 1, 1, 1],
                     [1, N, 1, 1],
                     [1, 1, N, 1],
                     [1, 1, 1, N]])*(mu/4))

k = 2 # Transition-transversion ratio
K80 = Model('K80') # Kimura 1980
K80.addRM(np.array([[N, k, 1, 1],
                    [k, N, 1, 1],
                    [1, 1, N, k],
                    [1, 1, k, N]]))

models = [JC69, K80]