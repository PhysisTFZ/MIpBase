# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 19:52:29 2022

@author: withp
"""
import csv
import numpy as np
def ReadConstrData(path, M, NK, initialCol):
    with open(path, encoding='gb18030', errors='ignore') as f:
        file = csv.reader(f)
        # each row of the table is a element of the data
        data = list(file)
        # number of the bases
        length = len(data)
        para = np.ones((length - 1, NK, M))
        curCol = col
        for nk in range(NK):
            for m in range(M):
                # skip the head of the table
                for i in range(1, length):
                    # assertion
                    if curCol >= len(data[i]) :
                        raise Exception('size error')
                    para[i - 1, nk, m] = data[i][curCol]
                # skip to the next col
                curCol += 1         
        return para, curCol

def ReadMinMaxPara(path, M, NK, initialCol):
    with open(path, encoding='gb18030', errors='ignore') as f:
        file = csv.reader(f)
        # each row of the table is a element of the data
        data = list(file)
        # number of the bases
        length = len(data)
        paraMin = np.ones((length - 1, NK, M))
        paraMax = np.ones((length - 1, NK, M))
        curCol = col
        for nk in range(NK):
            for m in range(M):
                # skip the head of the table
                def Loop (para, nk, m, curCol):
                    for i in range(1, length):
                        # assertion
                        if curCol >= len(data[i]) :
                            raise Exception('size error')
                            paraMax[i - 1, nk, m] = data[i][curCol]
                return (curCol + 1)
                curCol = Loop(paraMax, nk, m, curCol)
                curCol = Loop(paraMin, nk, m, curCol)
        return paraMax, paraMin, curCol

col = 1    
NK = 4
M = 3
path = 'data.csv'
col = 1
c1, col =  ReadConstrData(path,  M, NK, col)  
c2, col =  ReadConstrData(path,  M, NK, col)  
c3, col =  ReadConstrData(path,  1, 1, col)    
c4, col =  ReadConstrData(path,  M, NK, col)       
c5, col =  ReadConstrData(path,  1, 1, col)  
c6, col =  ReadConstrData(path,  1, 1, col)  
path = 'data.csv'  
v1, col =  ReadConstrData(path,  M, NK, col)  
v2, col =  ReadConstrData(path,  M, NK, col)  
v3, col =  ReadConstrData(path,  1, 1, col)    
v4, col =  ReadConstrData(path,  M, NK, col)       
v5, col =  ReadConstrData(path,  1, 1, col)  
v6, col =  ReadConstrData(path,  1, 1, col)  
pMax, PMin, col = ReadMinMaxPara(path, M, NK, col)
pMax, PMin, col = ReadMinMaxPara(path, M, NK, col)