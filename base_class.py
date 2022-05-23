# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 09:40:21 2022

@author: withp
"""
import numpy as np
import random
import math
import copy
DIM = 2
NK = 4
SQUARE = 2
LN2 = 0.30102999566398
EPI = 1/ (2*0.25)
# store all the optimization model parameter    
class Parameter :
    def __init__ (self, M, N, neighbor): 
        self.N = N
        self.M = M
        self.U = []
        self.uPara = np.zeros(N)
        self.uMin = []
        self.uMax = []
        self.pMax  = []
        # first dimension gamma, second dimension V
        self.constr = []
        self.linearPara = np.zeros((N, DIM, M, NK))
        self.neighbor = neighbor
        self.neighborNum = NK*np.ones(N)
        # first dimension c1, second dimension v1
        self.cv1 = 0.6*np.ones((N, DIM, NK, M))
        self.cv2 = 0.3*np.ones((N, DIM, NK, M))
        self.cv3 = 0.3*np.ones((N, DIM, 1, 1))
        self.cv4 = 0*np.ones((N, DIM, NK, M))
        self.cv5 = 0.05*np.ones((N, DIM, 1, 1))
        sign = [1, -1]
        for i in range(DIM) :
            for n in range(N):
                self.cv5[n,i,0,0] =  random.choice(sign) * self.cv5[n,i,0,0]
                
        self.cv6 = np.zeros((N, DIM, 1, 1))
    
    
    # get linear parameter of the model
    def GetLinearPara(self) :
        return self.linearPara
    
    # initialize the user parameter we used in the model
    def UserParaIntial(self, U, uMin, uMax):
        self.U = U
        self.uMin = uMin
        self.uMax = uMax
        for n in range(self.N):
            self.uPara[n] = (self.U[n] - self.uMin[n])/self.uMax[n] 
                
    def PowerParaIntial(self, pMax):
        self.pMax = pMax
        
    def ModelParaIntialC(self, c1, c2, c3, c4, c5, c6):
        self.cv1[0] = c1
        self.cv2[0] = c2
        self.cv3[0] = c3
        self.cv4[0] = c4
        self.cv5[0] = c5
        self.cv6[0] = c6
        
    def ModelParaIntialV(self, v1, v2, v3, v4, v5, v6):
        self.cv1[1] = v1
        self.cv2[1] = v2
        self.cv3[1] = v3
        self.cv4[1] = v4
        self.cv5[1] = v5
        self.cv6[1] = v6
    
    def SqaureConst(self, n, k, m, dim) :
        return (self.cv1[n, dim, k, m])
    def LinearConst(self, n, k, m, dim) :
        return (self.cv2[n, dim, k, m]*self.uPara[n])
    def ExpConst(self, n, k, m, dim) :
        e1 = -math.exp(-EPI *((1 + self.uPara[n])**SQUARE))
        slope = e1 + 1
        return (self.cv4[n, dim, k, m]*slope)     
    
    
    def ConstrParaInitial(self, constr):
        self.constr = copy.deepcopy(constr)
        # constr is a three dimension matrix
        for i in range(constr.shape[0]):
            for n in range(constr.shape[1]):
                    self.constr[i, n] -= self.cv3[i, n, 0, 0]*self.uPara[n]
                    self.constr[i, n] -= self.cv6[i, n, 0, 0]
    
    
    # Initialize the linear parameter of the constraints for the non-log form
    def LinearParaInitial(self) :
        # from the speed to volume
        for n in range(self.N):
            # from the speed to volume
            for i in range(DIM):
                for m in range(self.M):                    
                    # for each neighbor of the n
                    for k in range(int(self.neighborNum[n])) :
                        self.linearPara[n, i, m, k] = (
                                self.SqaureConst(n, k, m, i) + 
                                self.LinearConst(n, k, m, i) +
                                self.ExpConst(n, k, m, i))
                        
    # Initialize the constraints parameter with linear function
    def LinearConstrInitial(self, constr):
        self.constr = copy.deepcopy(constr)
        # constr is a three dimension matrix
        for n in range(constr.shape[0]):
            for i in range(constr.shape[1]):
                    self.constr[n, i] -= self.cv3[n, i, 0, 0]*self.uPara[n]
                    self.constr[n, i] -= self.cv6[n, i, 0, 0] 
        # for exp function, we need to plus the value at 0
                    for k in range(int(self.neighborNum[n])) :
                        for m in range(self.M):
                            self.constr[n, i] += self.cv4[n, i, k, m] 
     
    def Initial(self, pMax, U, uMin, uMax, constr) :
        self.UserParaIntial(U, uMin, uMax)
        self.PowerParaIntial(pMax)
        self.LinearConstrInitial(constr)
        self.LinearParaInitial()
                              
    """
    Generate the linear parameter of the constraints for the non-log form
    """ 
    def GenerateLinearPara(self) :
        # initial result is {1} for all the variables
        result = np.ones((self.N, self.M))
        # from the speed to volume
        for n in range(self.N):
            gradient = []
            residual = self.PositiveLinearLog(result, n, gradient)
            for i in range(DIM):
                self.constr[n, i] -= residual * self.cv6[n, i, 0, 0] 
                for m in range(self.M):                    
                    # for each neighbor of the n
                    for k in range(int(self.neighborNum[n])) :
                        self.linearPara[n, i, m, k] = (
                                self.SqaureConst(n, k, m, i) + 
                                self.LinearConst(n, k, m, i) +
                                self.ExpConst(n, k, m, i) + 
                                self.cv6[n, i, 0, 0] *gradient[m*k])
                        

                        
        
        