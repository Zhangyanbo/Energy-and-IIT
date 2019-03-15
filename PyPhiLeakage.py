from sympy.combinatorics.graycode import GrayCode
from itertools import permutations
import numpy as np
import warnings
import pyphi
import pandas as pd
import time

ntimes = 10
each_time_compute = 20

warnings.filterwarnings('ignore')

labels = ('A', 'B', 'C')
m = [[0, 0, 0],
     [0, 0, 1],
     [1, 0, 1],
     [1, 0, 0],
     [1, 1, 0],
     [1, 1, 1],
     [1, 1, 1],
     [1, 1, 0]]

def is_in(state, tpm):
    for astate in tpm:
        if state == astate.tolist():
            return True
    return False

def getphi(tpm):
    network = pyphi.Network(tpm, node_labels=labels)
    phis = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                state = (i, j, k)
                node_indices = (0, 1, 2)
                if is_in(list(state), tpm):
                    subsystem = pyphi.Subsystem(network, state, node_indices)
                    phis += [pyphi.compute.phi(subsystem)]
    return phis


# Show CPU Leakage

for i in range(ntimes):
    start = time.time()
    for j in range(each_time_compute):
        temp = getphi(np.array(m))
    end = time.time()
    print('epco', i, 'time cost:', end-start, 's')
