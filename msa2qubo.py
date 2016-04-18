#!/usr/bin/env python3

import os
import argparse
import math
import numpy
import copy
import sys
from Bio import SeqIO
from bvc import BVC


class Simulator:
	"""This class is still work in progress.  The idea is to be able to simulate the D-Wave machine for a very small number
	of binary variables: < 30; using a brute force approachCalculates E0 from a defined set of binary variables"""

	__index = 0
	__bv = []

	def __init__(self, bvc):
		self.data = []
		self.__bvc = bvc
		self.__bv = [0] * bvc.get_NbBV()

	def __x(self, k, j):
		val = 0
		offset = (k * self.__bvc.K() + j) * self.__bvc.m()
		for i in range(self.__bvc.m()):
			val += (self.__bv[offset + i] * 2 ** i)
		return val

	def __G(self, k, j):
		val = 0
		offset = self.__bvc.get_gVarOffset() + (k * self.__bvc.K() + j) * self.__bvc.p()
		for i in range(self.__bvc.p()):
			val += (self.__bv[offset + i] * 2 ** i)
		return val

	def __calcE0(self):
		"""Calculates E0 from a defined set of binary variables"""
		x = self.__bvc.x()
		sum_k = 0
		for k in range(self.__N):
			klen = len(self.__x[k])
			sum_j = 0
			for j in range(klen - 1):
				term_1 = self.__x(k, j + 1) - self.__x(k, j) - 1 - self.__G(k, j + 1)
				term_2 = self.__x(k, 0) - self.__G(k, 0)
				sum_j += (term_1 ** 2) + (term_2 ** 2)

			sum_k += sum_j

		return self.__bvc.l0() * sum_k

	def __calcE1(self):
		"""Calculates E1 from a defined set of binary variables"""
		x = self.__bvc.x()
		sum_k = 0.0
		for k in range(self.__bvc.N()):
			sum_j = 0.0
			for j in range(1, len(x[k])):
				sum_j += self.__G(k, j) ** 2
			sum_k += sum_j

		return self.__bvc.l1() * sum_k

	def __calcE2(self):
		"""Calculates E2 from a defined set of binary variables"""
		x = self.__bvc.x()
		N = self.__bvc.N()
		sum_k = 0.0
		for k in range(N - 1):
			sum_q = 0.0
			for q in range(k + 1, N):
				sum_i = 0.0
				for i in range(len(x[k])):
					sum_j = 0.0
					for j in range(len(x[q])):
						w_ijkq = self.__bvc.w(i,j,k,q)
						r_ijkq = self.__bvc.r(i,j,k,q)
						if w_ijkq != 0.0 and r_ijkq != 0.0:
							sum_j += w_ijkq * r_ijkq * math.floor(1 - self.__bvc.l2() * ((self.__x(k, i) - self.__x(q, j)) ** 2))
					sum_i += sum_j
				sum_q += sum_i
			sum_k += sum_q

		return -sum_k

	def __calcSolutionForInstance(self):
		"""Calculates Solution value from a particular instance of defined binary variables"""
		return self.__calcE0() + self.__calcE1() + self.__calcE2()

	def __incrementBV(self):
		"""Increment the state of binary variables by 1"""
		bvlen = len(self.__bv)
		bvstr = format(self.__index, '0' + str(bvlen) + 'b')
		for i in range(bvlen):
			self.__bv[i] = int(bvstr[i])

		# Value for next iteration
		self.__index += 1

	def calcSolution(self):
		"""Calculates optimal solution by trying all possible permutations of binary variables"""
		min_solution = copy.deepcopy(self.__bv)
		min_val = sys.float_info.max
		nbPermutations = 2 ** len(self.__bv)
		for i in range(nbPermutations):
			self.__incrementBV()
			val = self.__calcSolutionForInstance()
			if (val < min_val):
				min_val = val
				min_solution = copy.deepcopy(self.__bv)

		return min_val, min_solution


def main():
	parser = argparse.ArgumentParser("Convert Fasta file containing multiple sequences to align, into QUBO format")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output", required=True,
						help="The output file, containing the QUBO representation of the MSA problem")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA")
	parser.add_argument("-d", "--delta", type=float, default=2.0, help="Delta.  The scaling factor to apply when converting products of 3BVs to 2BVs.")
	parser.add_argument("-s", "--simulate", action='store_true', default=False,
						help="Whether to try and simulate D-Wave and come up with a solution to the problem.  Only runs if number of binary variables is < 30.")
	parser.add_argument("-l0", "--position_weighting", type=float, default=0.8,
						help="The weighting to apply to positioning of elements (must be larger than --gap_weighting)")
	parser.add_argument("-l1", "--gap_weighting", type=float, default=0.1,
						help="The weighting to apply to gap penalties")
	parser.add_argument("-l2", "--reward_weighting", type=float, default=10.0,
						help="The weighting to apply to reward matches (must be greater than 1.0)")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	if args.position_weighting < args.gap_weighting:
		print("-l0 must be larger than -l1")
		exit(1)

	if args.gap_weighting <= 0.0:
		print("-l1 must be larger than 0.0")
		print(1)

	if args.reward_weighting <= 1.0:
		print("-l2 must be > 1.0")
		exit(1)

	if args.delta <= 1.0:
		print("-d must be > 1.0")
		exit(1)

	print("Loading input into memory...", end="")
	handle = open(args.input, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	print(" done")
	print()
	print("Input:")

	# Collect variables
	bvc = BVC(P=args.P, d=args.delta, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting)
	for r in records:
		bvc.add_record(r)

		# Input will only be small so just output every sequence found
		print(">" + r.id)
		print(r.seq)

	print()

	# Print current state of variables
	print(bvc)

	# Save settings file
	print()
	bvc.save_settings(args.output + ".settings")
	print("Saved settings to " + args.output + ".settings")

	# Create matrix
	print()
	print("Creating W matrix ...", end="")
	sys.stdout.flush()
	bvc.createW()
	# print (m)
	print(" done")
	print("Creating matrix of coefficients of binary variables ...", end="")
	sys.stdout.flush()
	bvc.createBVMatrix()
	# print (m)
	print(" done")
	print("Number of active binary variables: ", bvc.calcActiveBVs())

	# Write QUBO file to disk
	print("Writing QUBO output to disk ...", end="")
	sys.stdout.flush()
	bvc.writeQUBO(args.output, args.input)
	print(" done")
	print("Output saved to: " + args.output)
	sys.stdout.flush()

	# Find solution???
	if args.simulate:
		if bvc.get_NbBV() > 30:
			print("Not trying to simulate data... way too many binary variables to finish in a reasonable time!",
				  file=sys.stderr)
		else:
			sim = Simulator(bvc)
			val, bv = sim.calcSolution()

			print("Min value is: ", val)
			print("Binary variables are:\n", bv)


main()
