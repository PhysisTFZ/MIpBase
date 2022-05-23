# -*- coding: utf-8 -*-
"""
Created on Sun May 15 19:35:07 2022

@author: withp
"""
import base_class as base
from pyscipopt import Model, quicksum, log, exp    
import time         
class SolverModel(base.Parameter) :
    def __init__ (self, M, N, neighbor): 
        super(SolverModel, self).__init__(M, N, neighbor)
        self.M = M
        self.N = N
        # neighbor is two dimension array
        self.neighbor = neighbor
        self.model = Model("wireless")
        self.bp = {}
        self.result = base.np.zeros((N, M))
       
    def GetResult(self) :
        return self.result
    
    # build the variables   
    def ModelIntial(self):
        for n in range(self.N):
            for m in range(self.M):
                self.bp[n, m] = self.model.addVar(vtype="B", name="bp/%s/%s"%(n, m)) 
   
    # generate the objective function 
    def GenObjective(self):
        self.model.setObjective(quicksum(self.bp[n, m]*self.pMax[n, m]  for (n, m) in self.bp), 'minimize')  
   
    def SqaureConstraits(self, n, k, m, dim) :
        id = int(self.neighbor[n, k])
        return (self.cv1[n, dim, k, m]*(self.bp[id, m])**base.SQUARE)
    def LinearConstraits(self, n, k, m, dim) :
        id = int(self.neighbor[n, k])
        return (self.cv2[n, dim, k, m]*self.uPara[n]*self.bp[id, m])
    def ExpConstraits(self, n, k, m, dim) :
        id = int(self.neighbor[n, k])
        return (self.cv4[n, dim, k, m]* exp(-base.EPI *(self.bp[id, m] *((1 + self.uPara[n])**base.SQUARE)))) 
    # special case, the sum is inside the log function
    def LogConstraits(self, n, dim) :
        constr = 0
        ratio = 1/(self.M*(1 +  self.neighborNum[n] ))
        for m in range(self.M):
            for k in range(base.NK):
                id = int(self.neighbor[n, k])
                constr += ratio*exp(1/2 * self.bp[id, m] *  (self.uPara[n] + 1 ))
        return log(constr)       
    
    def GenConstraints(self) :
        for i in range(base.DIM):
            for n in range(self.N):
                    constr = 0
                    for m in range(self.M):
                        for k in range(base.NK):
                            id = int(self.neighbor[n, k])                        
                            constr +=  self.bp[id, m] * self.linearPara[n, i, m, k]
                    constr += self.LogConstraits(n, i)            
                    self.model.addCons(constr >= self.constr[n, i])
                    
    def ModelBuild(self):
        self.ModelIntial()
        self.GenConstraints()
        self.GenObjective()
        
    def SolveModel(self):
        self.model.getObjectiveSense()
        start = time.time()
        self.model.optimize()
        elapsed = (time.time() - start)
        print("Scip time used: " + str(elapsed))
        stage = self.model.getStatus()
        if stage == 'optimal':
            cost = self.model.getBestSol()
            print("The optimal value is", self.model.getObjVal()) 
            for v in self.model.getVars():
                n =  int(v.name.split('/')[1])
                m =  int(v.name.split('/')[2])
                value = self.model.getVal(v)
                if value > 0.001:
                    self.result[n, m] = value
                    #print(v.name, "=", value)
            return cost
        else: 
            return 0
    
# optimization model class, inherit from the Paramater class        
class OptiModel(SolverModel): 
        # model building, with M, N, T, neightbor information
        def __init__ (self, M, N, neighbor): 
                super(OptiModel, self).__init__(M, N, neighbor)
                
        def ParaIntial (self, U, uMin, uMax, pMax, constrPara):
            super(OptiModel, self).Initial(pMax, U, uMin, uMax, constrPara)

            
        def ModelInitial(self):
            super(OptiModel, self).ModelBuild()

