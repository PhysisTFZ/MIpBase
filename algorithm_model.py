# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:09:45 2022

@author: withp
"""

from constraint import Constraint
from constraint import NonZero
import random
import base_class as base
import cvxpy as cp
import time
class AlgorithmModel:
    def __init__ (self,N,M) :
        self.M = M
        self.N = N
        self.neighbor = []
        self.pMax = []
        self.constraints = []
        self.maxIter  = 3
        self.iter = 1
        self.optimalResult = []
    
    def ResultStore(self, result):
        for n in range(len(self.constraints)) :
            self.constraints[n].AppendResult(result)
     
    # initial the constraints
    def InitialConstraints(self, result):
        for n in range(self.N):
         self.constraints[n].Initial()
         self.constraints[n].AppendResult(result)
         
    # initialize the total result            
    def Initial(self, parameter):
        initial = base.np.zeros((self.N, self.M))
        self.pMax = parameter.pMax 
        self.neighbor = parameter.neighbor
        for n in range(self.N):
            self.constraints.append(Constraint(n, parameter))
            self.constraints[n].AppendResult(initial)
     
            
    def ConstraintsBuild(self, constraints, variable) :        
        for n in range(len(self.constraints)) :
            for i in range(base.DIM) : 
                linear = []
                value = []
                neighborNum  = self.constraints[n].GetNeighborNum()
                self.constraints[n].GenearteLinearPara(i, linear, value)
                for l in range(len(linear)) :
                    constr = 0
                    for m in range(self.M) :
                        for k in range(neighborNum):
                            id = int(self.neighbor[n, k])                        
                            constr +=  (variable[id, m] * linear[l][m, k])
                    constraints += [constr >=  value[l]]
            
                    
    def ModelBuild(self) :
        constraints = []
        b = cp.Variable((self.N, self.M), boolean=True)
        self.ConstraintsBuild(constraints, b)
        cost = 0
        for n in range(self.N):
            for m in range(self.M):
                cost +=  b[n, m] * self.pMax[n, m] 
        objects = cp.Minimize(cost)
        self.model = cp.Problem(objects, constraints)       
        start = time.time()    
        self.model.solve('SCIP')
        #print("constraints length is :" + str(len(constraints)))
        elapsed = (time.time() - start)
        #print("Cvxpy time used: " + str(elapsed))
        #print("The optimal value is", self.model.value) 
        self.ResultStore(b.value)
        #print(b.value)
        self.optimalResult = b.value
        return self.model.status
        
    def ConstraintsCheck(self, resultId):
        legal = True
        for n in range(len(self.constraints)) :
            for i in range(base.DIM) :   
                result = self.constraints[n].ConstraintsCheck(resultId, i)
                if result[0] == False:
                    """print(" model value: " + str(result[1]) +" real value: " +
                          str(result[2]))"""
                    legal = False
                    self.constraints[n].AddExpansionId(i, resultId)
        return legal
    
     
    def Iteration(self):
        for iter in range(self.maxIter) :
            if self.ConstraintsCheck(iter):
                cost = 0
                for n in range(self.N):
                    for m in range(self.M):
                        cost +=  self.optimalResult[n, m] * self.pMax[n, m] 
                print("Optimal value is: " + str(cost) +" Iteration is : " + str(iter)) 
                return True
            if self.ModelBuild() == "infeasible" :
                print("Infeasible and Iteration is : " + str(iter))
                return True
        return False
    
    def Improve(self, target) :
        for m in range(self.M) :
            for n in range(self.N) :
                neighborNum  = self.constraints[n].GetNeighborNum()
                one = []
                for k in range(neighborNum):
                    id = int(self.neighbor[n, k])  
                    if NonZero(self.optimalResult[id, m]) == False :
                        one += [k] 

                if (neighborNum - len(one)) < target :
                    choice = target + len(one) - neighborNum
                    result = random.sample(one, choice)
                    for i in range(len(result)) :
                        id = int(self.neighbor[n, result[i]])     
                        self.optimalResult[id, m] = 1  
                    
                    
                    
    def Process(self) :
        while self.Iteration() == False :
            self.Improve(self.iter)
            #resulr = 0.5*(one + self.optimalResult)
            self.InitialConstraints(self.optimalResult) 
            self.iter += 1
            
    def Algorithm(self, parameter) :  
        self.Initial(parameter)
        start = time.time()    
        self.Process()
        elapsed = (time.time() - start)
        print("Algorithm time used: " + str(elapsed))
        
        
           