"""
    created by tfz 2022.3.15, used for Mip test about the article model
    of "Deep Reinforcement Learning for Joint Channel
    Selection and Power Control in D2D Networks"
"""
import gurobipy as grb
import base_class as base
import time
import math

def ModelBuild(para):   
    model = grb.Model("wireless_trans")
    b = model.addVars(para.N, para.M, vtype = grb.GRB.BINARY, name = "b")  
    x = model.addVars(para.N, vtype = grb.GRB.CONTINUOUS, name = "x")  
    y = model.addVars(para.N, vtype = grb.GRB.CONTINUOUS, name = "y")  
    model.update()
    
    for n in range(para.N):
        ratio = 1/(para.M *(1 + base.NK))
        e1 = math.exp(1/2 *((1 + para.uPara[n])))
        slope = ratio* (e1 - 1)
        bSum = 0
        for m in range(para.M):
            for k in range(base.NK):
                id = int(para.neighbor[n, k])                        
                bSum += slope* b[id, m] +  ratio
        model.addConstr (x[n] == bSum)     
        model.addGenConstrExp(y[n], x[n])         
        for i in range(base.DIM):
            constr = 0
            for m in range(para.M):
                for k in range(base.NK):
                    id = int(para.neighbor[n, k])                        
                    constr += b[id, m] * para.linearPara[n, i, m, k]
            constr += para.cv5[n,i,0,0] *y[n] 
        model.addConstr (constr >= para.constr[n, i])       
    model.write("ss.lp")
    model.setObjective(grb.quicksum(para.pMax[n, m] * b[n, m] for (n, m) in b), grb.GRB.MINIMIZE)
    model.setParam("NonConvex", 2)
    return model
  
def ModelSolve(model):         
  model.setParam('OutputFlag', 0)
  start = time.time()
  model.optimize()
  elapsed = (time.time() - start)
  print("Gurobi time used: " + str(elapsed))
  if (model.status == grb.GRB.OPTIMAL):
      print("Gurobi Model Cost %s"%model.objVal)
  else :
      print("failure")
  """for i in model.getVars():
      print('%s = %g' % (i.varName, i.x), end = " ")
  return"""