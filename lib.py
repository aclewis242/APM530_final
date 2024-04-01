### LIBRARY OF MODELS

from model import *
import numpy as np

JC69 = Model('JC69') # Jukes & Cantor 1969 -- simplest model, assumes equal mutation rates
JC69.rate_mat = np.array([[-3, 1, 1, 1],
                          [1, -3, 1, 1],
                          [1, 1, -3, 1],
                          [1, 1, 1, -3]])*(JC69.sub_rate/4)