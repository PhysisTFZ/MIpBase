# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 20:10:44 2022

@author: withp
"""
import numpy as np
import random 
DIM = 2
VMAX = 10
GAMMAMAX = 10
MUL = 0.05
PMAX = 10
UMIN = 30
UMAX = 60
thetaMax = 9
def GenNeighbor(N, neighNum):
    neighbor = np.ones((N, neighNum))
    for i in range(N):
        result_list = list(range(0, N))
        result_list.remove(i)
        np.random.shuffle(result_list)
        neighbor[i, 0]  = i
        for j in range(1, neighNum):
            neighbor[i, j] =  result_list[j - 1]
    return neighbor
    
def GenUPara(M, N) :
    U = np.ones(N)
    uMin = np.ones(N)
    uMax = np.ones(N)
    for n in range(N):
        uMax[n] = random.randrange(UMIN + 1, UMAX, 1) 
        uMin[n] = random.randrange(1, UMIN, 1) 
        U[n] = random.randrange(uMin[n] + 1 , uMax[n] + 1, 1) 
    return U, uMin, uMax

def GenPPara(M, N) :
    pMax = np.ones((N, M))
    for n in range(N):
        maxV =  random.randrange(1, PMAX, 1) 
        for m in range(M):
            pMax[n, m] = 1#maxV
    return pMax

def GenConPara(N) :
    constr = np.ones((N, DIM))
    for n in range(N):
        for i in range(DIM):
            constr[n, i] = 0.1 * random.randrange(1, GAMMAMAX, 1) 
            constr[n, i] = 0.1 * random.randrange(1, VMAX, 1) 
    return constr
