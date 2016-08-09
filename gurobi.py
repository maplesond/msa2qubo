#!/usr/bin/env python3

import bvc

from gurobipy import *
import numpy as np


def optimise(data):

	# Create a new model
	m = Model("qp")

	# Create variables
	x_k = 0
	G_k = data.get_gVarOffset(intmode=True)
	vars = [None] * data.get_NbIV()
	for k in range(data.N()):
		L_k = data.lenK(k)
		for j in range(L_k):
			x_kj = x_k + j
			vars[x_kj] = m.addVar(name="x_" + str(k) + "," + str(j), vtype=GRB.INTEGER)
		x_k += L_k

	for k in range(data.N()):
		L_k = data.lenK(k)
		for j in range(L_k):
			G_kj = G_k + j
			vars[G_kj] = m.addVar(name="G_" + str(k) + "," + str(j), vtype=GRB.INTEGER)
		G_k += L_k


	# Integrate new variables
	m.update()

	data.createBVMatrix(intmode=True)

	print("Quadratic coefficients:")
	data.printqim()
	print()
	print("Linear coefficients:")
	data.printlil()
	print()

	# Set objective: x^2 + x*y + y^2 + y*z + z^2 + 2 x
	obj = data.ienergy
	for i in range(data.get_NbIV()):
		for j in range(data.get_NbIV()):
			if data.qim(i,j) != 0:
				obj += data.qim(i,j) * vars[i] * vars[j]
		if data.lil(i) != 0:
			obj += data.lil(i) * vars[i]

	obj = data.l0() * (obj)
	print("Integer Objective Function:")
	print(obj)
	print()
	m.setObjective(obj)


	for i in range(data.get_NbPositioningVars(intmode=True)):
		m.addConstr(vars[i] >= 0, "cx" + str(i))
		m.addConstr(vars[i] <= data.M(), "cx" + str(i))

	for i in range(data.get_gVarOffset(intmode=True),data.get_NbIV()):
		m.addConstr(vars[i] >= 0, "cG" + str(i-data.get_gVarOffset(intmode=True)))
		m.addConstr(vars[i] <= data.P(), "cG" + str(i-data.get_gVarOffset(intmode=True)))


	m.optimize()

	for v in m.getVars():
		print('%s: %g' % (v.varName, v.x))

	print('Obj: %g' % obj.getValue())
