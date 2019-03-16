from sympy.combinatorics.graycode import GrayCode
from itertools import permutations
import numpy as np
import warnings
import pyphi
import pandas as pd
import time
import sys

start = time.time()

################
# Terminal input:
# **.py start end total filename
# Note: the end is not included
################

start_index = int(sys.argv[1])
end_index = int(sys.argv[2])
total_num = int(sys.argv[3])
phi_path = './data/' + sys.argv[4]

print('Compute Phi from label', start_index, 'to', end_index)
print(end_index - start_index, 'compute of', total_num, 'all computes.')

#total_compute = 100
astep = (end_index - start_index) // 10

ising_path = './data/isingEnergy_0_10080.csv'
#phi_path = './data/Phi_0_10080.csv'

##############Functions##############
def ReverseBinCode(l):
    # Generate reversed binary code for pyphi
    maxIndex = 2 ** l - 1
    return [np.binary_repr(maxIndex - index, width=l) for index in range(2 ** l)]



def code2str(codeList):
    s = ''
    for i in codeList:
        s += str(i)
    return s

def relabel(m0, order):
    length = len(m0[0])
    bcode = ReverseBinCode(length)
    reorderRule = {}
    for originLabel, newLabel in zip(bcode, order):
        reorderRule[originLabel] = newLabel
    #print(reorderRule)
    rule = {}
    # 初始化rule
    for currentState, nextState in zip(bcode, m0):
        rule[code2str(currentState)] = nextState
    # 更新顺序
    for currentState, nextState in zip(bcode, m0):
        #print(reorderRule[code2str(currentState)], 'to', str(nextState))
        rule[code2str(reorderRule[code2str(currentState)])] = reorderRule[code2str(nextState)]
    
    new_m = [rule[startState] for startState in ReverseBinCode(length)]
    
    return np.array(new_m)

def AllLabelConfigurations(mat):
    # 返回所有可能的重标记方法
    _check = {}
    all_elements = [
    [1,1,1],
    [1,1,0],
    [1,0,1],
    [1,0,0],
    [0,1,1],
    [0,1,0],
    [0,0,1],
    [0,0,0]
    ]
    for order in permutations(all_elements):
        _relabel = relabel(mat, list(order))
        _check[str(_relabel)] = _relabel
    return list(_check.values())

# 输入状态转移列表，给出各个量的卡诺图矩阵

def MapToKarnaughMap(_map):
    # Turn map to MapToKarnaugh Maps
    # Generate ReverseBinCode
    index = ReverseBinCode(len(_map[0]))
    #index = np.array([[int(char) for char in strcode] for strcode in index])
    
    dim1 = (int(len(_map[0]) / 2))
    dim2 = (len(_map[0]) - int(len(_map[0]) / 2))
    #print('dim =', dim1, 'x', dim2)
    
    gc1 = GrayCode(dim1)
    gc2 = GrayCode(dim2)
    gd1 = dict(zip(list(gc1.generate_gray()), range(2 ** dim1)))
    gd2 = dict(zip(list(gc2.generate_gray()), range(2 ** dim2)))
    #print(gd1)
    #print(gd2)
    
    _map = np.array(_map).T
    KarnaughMaps = []
    for targetStates in _map:
        #print('-------')
        Karnaugh = np.array(([[0] * (2 ** dim1)]) * (2 ** dim2)).T
        #print(Karnaugh)
        for state, mcode in zip(index, targetStates):
            #print(state, '->', mcode)
            #print('c1=',state[0:dim1])
            #print('c2=',state[dim1:])
            Karnaugh[gd1[state[0:dim1]], gd2[state[dim1:]]] = mcode
        KarnaughMaps += [Karnaugh]
        #print(Karnaugh)
    
    return np.array(KarnaughMaps)

def IsingEnergy(mat, energy={0:-1, 1:1}):
    # 给出mat的Ising能量
    # Return Ising energy of mat
    mat = np.array(mat)
    # If mat is a list of mat, then return all energy
    if len(mat.shape) == 3:
        return np.array([IsingEnergy(_mat, energy) for _mat in mat])
    # Otherwise return this energy
    x, y = mat.shape
    total_E = 0
    for i in range(x):
        for j in range(y):
            total_E += energy[mat[i][j]] * energy[mat[(i - 1) % x][j]]
            total_E += energy[mat[i][j]] * energy[mat[(i + 1) % x][j]]
            total_E += energy[mat[i][j]] * energy[mat[i][(j - 1) % y]]
            total_E += energy[mat[i][j]] * energy[mat[i][(j + 1) % y]]
    return -total_E

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


###################

m = [[0, 0, 0],
     [0, 0, 1],
     [1, 0, 1],
     [1, 0, 0],
     [1, 1, 0],
     [1, 1, 1],
     [1, 1, 1],
     [1, 1, 0]]

all_labels = AllLabelConfigurations(m)
print('total labels:',len(all_labels))

warnings.filterwarnings('ignore')

labels = ('A', 'B', 'C')

def IsingEnergyExp(exp_m, show_process=True, step = 100, saveQ=False):
    print('Computing Ising Energy...')
    #ising_energy = [IsingEnergy(MapToKarnaughMap(alabel)) for alabel in exp_m]
    ising_energy = []
    for st in range(0, len(exp_m), step):
        ising_energy += [IsingEnergy(MapToKarnaughMap(alabel)) for alabel in exp_m[st:st+step]]
        #print('ising energy:', st / len(exp_m) * 100, '%')
    print('done.')
    if saveQ:
        isingEnergyData = pd.DataFrame(ising_energy)
        isingEnergyData.to_csv(ising_path)
    
    print('Computing Phi...')
    #phis = [getphi(alabel) for alabel in exp_m]
    phis = []
    for st in range(0, len(exp_m), step):
        phis += [getphi(alabel) for alabel in exp_m[st:st+step]]
        print('Phi:', (st + start_index + step) / total_num * 100, '%,', (st + step) / len(exp_m) * 100, '%','of this part')
        if saveQ:
            PhiData = pd.DataFrame(phis)
            PhiData.to_csv(phi_path)
    print('done.')
    
    return ising_energy, phis

def allEnergy(exp_m):
    temp = [IsingEnergy(MapToKarnaughMap(alabel)) for alabel in exp_m]
    return temp

#isingEnergys, phis = IsingEnergyExp(all_labels[start_index : end_index], show_process=True, step=astep, saveQ=True)

#PhiData = pd.DataFrame(phis)
#PhiData.to_csv(phi_path)

isingEnergys = allEnergy(all_labels)
isingEnergyData = pd.DataFrame(isingEnergys)
isingEnergyData.to_csv('./data/isingEnergy_0_10080.csv')

end = time.time()

print('total time cost:', end-start, 's')
