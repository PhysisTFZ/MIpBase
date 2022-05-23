# -*- coding: utf-8 -*-
"""
Created on Thu May 19 17:50:56 2022

@author: withp
"""
import numpy as np
import math
ACCURACY = 0.001

def NonZero(num) :
    if num < -ACCURACY:
        return True
    elif num > ACCURACY:
        return True
    else:
        return False
    
class Constraint :
    def __init__ (self, index, paramter) :
        self.index = index
        self.M = paramter.M
        self.uPara = paramter.uPara[index]
        # two form constraints, first : volume, second: speed
        self.linearPara = paramter.linearPara[index]
        # neighbor num of the index 
        self.neighborNum = paramter.neighborNum[index]
        # neighbor array of the index
        self.neighbor = paramter.neighbor[index]
        self.resultLength  = 0
        self.expansionId = [[], []]
        self.results = []
        self.constraint = paramter.constr[index]
        self.cv1 = paramter.cv1[index]
        self.cv2 = paramter.cv2[index]
        self.cv3 = paramter.cv3[index]
        self.cv4 = paramter.cv4[index]
        self.cv5 = paramter.cv5[index]
        self.cv6 = paramter.cv6[index]
    
    def Initial(self) :
        self.resultLength  = 0
        self.expansionId = [[], []]
        self.results = []
        
    # get the expansion Id
    def GetResultLength(self) :
        return int(self.resultLength)
    
    #get the cv5 value of the given dim
    def GetCV5(self, dim) :
        return int(self.cv5[dim, 0, 0] )
    #
    def AddExpansionId(self, i, id) :
        self.expansionId[i] += [id]
    # Increase the expansionId
    def IncreaseResultLength(self) :
        self.resultLength += 1
    # get neighbor num of 
    def GetNeighborNum(self) :
        return int(self.neighborNum)
    # append the new result   
    def AppendResult(self, result):
        self.results +=[result]
    # get inner sum of the log function in the constratints
    def GetInnerValueOfLog(self, result):
        ratio = 1/(self.M*(1 + self.neighborNum))
        value = 0
        for m in range(self.M):
            for k in range(int(self.neighborNum)):
                id = int(self.neighbor[k])
                value += math.exp(1/2 * result[id, m] * (self.uPara + 1))
        return (ratio*value) 
    
    # get inner sum of the log function in the constratints only consider neighbor
    def GetNeighborSum(self, neighbor):
        ratio = 1/(self.M*(1 + self.neighborNum))
        value = 0
        for m in range(self.M):
            for k in range(int(self.neighborNum)):
                value += math.exp(1/2 * neighbor[m, k] * (self.uPara + 1))
        return (ratio*value)                 
          
    def PositiveLinearLog(self, result, gradient):
        residual = 0
        value = self.GetInnerValueOfLog(result)
        logValue = math.log(self.GetInnerValueOfLog(result))
        reciprocal = 1/value
        ratio = 1/(self.M*(1 + self.neighborNum))
        for m in range(self.M) :
            for k in range(int(self.neighborNum)):
                id = int(self.neighbor[k])
                partialDiff = reciprocal* ratio*(math.exp(1/2 * result[id, m] *
                              (self.uPara + 1)) * 1/2 *(self.uPara + 1))
                gradient += [partialDiff]
                residual -= partialDiff*result[id, m]
        return (residual + logValue)  
    
    # find an interpolated point inside
    def NegtiveLinearLog(self, i, gradient, value) :
        first = np.zeros((self.M, int(self.neighborNum)))
        last = np.ones((self.M, int(self.neighborNum)))
        innerSum = []
        innerSum+=[self.GetNeighborSum(first)]
        if len(self.expansionId[i]):
            current = self.expansionId[i][-1]
            innerSum +=[self.GetInnerValueOfLog(self.results[current])]
        innerSum+= [self.GetNeighborSum(last)]
        innerSum.sort()
        ratio = 1/(self.M*(1 + self.neighborNum))
        # log(x0) + slope*(x1 - x0)
        for i in range(1, len(innerSum)) :
            dist = innerSum[i] - innerSum[i-1]
            if NonZero(dist) :
                slope = math.log(innerSum[i]) - math.log(innerSum[i-1]) /dist
                # 1 + slope*(x)
                e1 = ratio * math.exp(1/2 *(1 + self.uPara))
                value += [math.log(innerSum[i-1]) + slope * ratio]
                # outer slope * inner slope
                slope *= (e1 - 1) * ratio
                gradient +=[slope]     
    
    # to eliminate the infeasible results
    def GetredundantCut(self, i, cut, oneNum) :
        nNum = int(self.neighborNum)
        for l in range(len(self.expansionId[i])) :
            count = 0
            rId = self.expansionId[i][l]
            cut += [np.zeros((self.M, nNum))]
            for m in range(self.M) :
                for k in range(nNum) :
                    id = int(self.neighbor[k])
                    if NonZero(self.results[rId][id, m]) :
                        cut[l][m, k] = -1
                        count += 1
                    else:
                        cut[l][m, k] = +1
            
            oneNum +=[1 - count]
    """
    Generate the linear parameter of the constraints for the  - log form
    """ 
    def LinearParaWithNegativeLog(self, i, linear, value) :
        gradient =[]
        nNum = int(self.neighborNum)
        self.NegtiveLinearLog(i, gradient, value)
        for j in range(len(gradient)) :
            linear += [gradient[j] * np.ones((self.M, nNum))]
            for m in range(self.M) :
                for k in range(nNum) :
                    linear[j][m, k] = (self.linearPara[i, m, k] 
                                       + self.cv5[i, 0, 0]*linear[j][m, k])
            value += [self.constraint[i] - value[j] * self.cv5[i, 0, 0]] 
        self.GetredundantCut(i, linear, value)
                    
    """
    Generate the linear parameter of the constraints for the  + log form
    """ 
    def LinearParaWithPositiveLog(self, i, linear, value) :
        for l in range(len(self.expansionId[i])) :
            nNum = int(self.neighborNum)
            linear +=[np.zeros((self.M, nNum))]
            # from the speed to volume
            gradient = []
            resultId = self.expansionId[i][l]
            residual = self.PositiveLinearLog(self.results[resultId], gradient)
            for m in range(self.M):                    
                # for each neighbor of the n
                for k in range(nNum) :
                    linear[l][m, k] = (self.linearPara[i, m, k] + 
                                    self.cv5[i, 0, 0] *gradient[m*nNum +k])
            value += [self.constraint[i] - residual * self.cv5[i, 0, 0]] 
                    
    """
    Generate the linear parameter of the constraints for the  Non log form
    """ 
    def LinearParaWithoutLog(self, i, linear, value) :
        nNum = int(self.neighborNum)
        linear +=[np.zeros((self.M, nNum))]
        for m in range(self.M) :
            for k in range(nNum) :
                linear[0][m, k] = (self.linearPara[i, m, k])
        value += [self.constraint[i]]   
        
    """
    Generate the linear parameter of the constraints
    """ 
    def GenearteLinearPara(self, i, linear, value) :
        if NonZero(self.cv5[i, 0, 0]) :
            if self.cv5[i, 0, 0] > 0:
                self.LinearParaWithPositiveLog(i, linear, value)
            else:
                self.LinearParaWithNegativeLog(i, linear, value)
        else:
                self.LinearParaWithoutLog(i, linear, value)
    
    # get constraint value of the nth function                 
    def ConstraintsValue(self, rId, dim):
        value = 0
        for m in range(self.M):                    
            # for each neighbor of the n
            for k in range(int(self.neighborNum)) :
                # get the k neighbor index
                id = int(self.neighbor[k])
                value += (self.linearPara[dim, m, k] * 
                          self.results[rId][id, m]) 
        innerValue = self.GetInnerValueOfLog(self.results[rId])
        # get value of log function 
        value += self.cv5[dim, 0, 0] * math.log(innerValue)
        return value 
    
    # list all the violated result
    def ConstraintsCheck(self, rId, dim) :
        value = self.ConstraintsValue(rId, dim)
        if value < self.constraint[dim] :
            return [False, value, self.constraint[dim]]
        else:
            return [True, value, self.constraint[dim]]
                            
