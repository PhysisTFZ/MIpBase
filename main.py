# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 20:08:19 2022

@author: withp
"""
import  scipy_model as Smodel
import  gurobi_model as Gmodel
import algorithm_model as algo
import parameter as para
from base_class import Parameter

M = 3
N = 255
neighNum = 4
neighbor = para.GenNeighbor(N, neighNum)
U, uMin, uMax = para.GenUPara(M, N)
pMax = para.GenPPara(M, N)
constrPara = para.GenConPara(N)
optmodel = Smodel.OptiModel(M, N, neighbor)
optmodel.ParaIntial(U, uMin, uMax, pMax, constrPara)
cost = optmodel.SolveModel()
gmodel = Gmodel.ModelBuild(optmodel)
Gmodel.ModelSolve(gmodel)
pMax = para.GenPPara(M, N)
constrPara = para.GenConPara(N)

parameter = Parameter(M, N, neighbor)
parameter.Initial(pMax, U, uMin, uMax, constrPara)
model = algo.AlgorithmModel(N, M)
model.Algorithm(parameter)