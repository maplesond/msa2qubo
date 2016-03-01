#!/usr/bin/env python3

import os
import argparse
import math
import numpy
import copy
import sys
from Bio import SeqIO


class BVC:
	"""This class represents the binary variables required for a given number of sequences
	To use the class initialise, and optionally override default parameters representing the max gap size and weightings,
	then add BioPython sequence records that will represent the MSA.  Finally call createBVMatrix to create the symettric
	matrix containing binary variable coefficients."""
	__N = 0
	__Lmax = 0
	__K = 0

	__x = []

	__bvm = numpy.zeros((0, 0))

	def __init__(self, P=1, l0=0.8, l1=1.0, l2=10.0):
		self.data = []
		self.__P = P
		self.__l0 = l0
		self.__l1 = l1
		self.__l2 = l2

	def __str__(self):
		return "N=" + str(self.__N) + "\tNumber of sequences\n" \
									  "Lmax=" + str(self.__Lmax) + "\tLength of longest sequence\n" \
																   "K=" + str(
			self.__K) + "\tTotal number of characters\n" \
						"M=" + str(self.M()) + "\tLength of each row in solution space\n" \
											   "m=" + str(
			self.m()) + "\tNumber of binary variables required to represent each position in solution space\n" \
						"P=" + str(self.__P) + "\tMaximum gap size\n" \
											   "p=" + str(
			self.p()) + "\tNumber of binary variables required to represent each gap\n\n" \
						"l0=" + str(self.__l0) + "\tPosition weighting\n" \
												 "l1=" + str(self.__l1) + "\tGap weighting\n" \
																		  "l2=" + str(
			self.__l2) + "\tReward weighting\n\n" \
						 "Solution space will contain " + str(self.calc_solutionSpaceSize()) + " cells\n\n" \
																							   "# Binary Variables required: " + str(
			self.calc_minBV()) + "-" + str(self.calc_maxBV()) + "\n\n"

	def add_record(self, record):
		"""Adds a biopython record and updates the state of this object accordingly"""
		l = len(record.seq)
		self.__Lmax = max(self.__Lmax, l)
		self.__x.append([0] * l)  # Add a list the size of the this sequence to the 'x' list
		self.__K += l
		self.__N += 1

	def N(self):
		return self.__N

	def K(self):
		return self.__K

	def P(self):
		return self.__P

	def M(self):
		"""Returns the length of each row in solution space"""
		return 2 * self.__Lmax + 1

	def m(self):
		"""Returns the number of binary variables required to represent each row in solution space"""
		return math.ceil(numpy.log2(self.M()))

	def p(self):
		"""Returns the number of binary variables required to represent each gap"""
		return math.ceil(numpy.log2(self.__P + 1))

	def x(self):
		return self.__x

	def l0(self):
		return self.__l0

	def l1(self):
		return self.__l1

	def l2(self):
		return self.__l2

	def calc_solutionSpaceSize(self):
		"""Size of the solution space"""
		return self.__N * self.M()

	def calc_minBV(self):
		"""Calculate lower bound on number of binary variables required"""
		return self.__K * (self.m() + self.p())

	def calc_maxBV(self):
		"""Calculate upper bound on number of binary variables required"""
		return self.calc_minBV() + math.ceil(self.__N * (self.__N - 1) * (self.__Lmax ** 2) * (1 + 2 * self.m()) / 2)

	def get_NbBV(self):
		"""Return the actual number of binary variables required for this problem"""
		# TODO Need to work out how to get the number of bvs required for E2, for now just assume maxBV
		return self.calc_maxBV()

	def get_gVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a gap"""
		return self.__K * self.m()

	def get_rVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_gVarOffset() + (self.__K * self.p())

	def __addE0Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E0"""
		x_bits = self.m()
		g_bits = self.p()
		g_offset = self.get_gVarOffset()
		for k in range(0, len(self.__x) - 1):
			sum_j = 0
			for j in range(0, len(self.__x[k])):
				x1_pos = k * (j + 1) * x_bits
				x0_pos = k * j * x_bits
				g_pos = g_offset + (k * (j + 1) * g_bits)
				for i in range(0, x_bits - 1):
					# TODO Not sure what to do with linear coefficients
					x0i = x0_pos + i
					x1i = x1_pos + i
					gi = g_offset + g_pos + i
					self.__bvm[x1i, x1i] += self.__l0
					self.__bvm[x1i, x0i] -= self.__l0
					# self.__bvm[x1i,?] -= self.__l0   # Linear
					self.__bvm[x1i, gi] -= self.__l0
					self.__bvm[x0i, x1i] -= self.__l0
					self.__bvm[x0i, x0i] += self.__l0
					# self.__bvm[x0i,?] -= self.__l0   # Linear
					self.__bvm[x0i, gi] += self.__l0
					self.__bvm[gi, x1i] -= self.__l0
					self.__bvm[gi, x0i] += self.__l0
					# self.__bvm[gi,?] -= self.__l0    # Linear
					self.__bvm[gi, gi] += self.__l0

				# Coefficients for second term (gap prior to first character)
				x_pos = k * 1 * x_bits
				g_pos = g_offset + (k * 1 * g_bits)
				for i in range(0, x_bits - 1):
					xi = x_pos + i
					gi = g_offset + g_pos + i
					self.__bvm[xi, xi] += self.__l0
					self.__bvm[gi, gi] += self.__l0
					self.__bvm[xi, gi] -= self.__l0
					self.__bvm[gi, xi] -= self.__l0

	def __addE1Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E1"""
		bits = self.p()
		offset = self.get_gVarOffset()  # The offset required to get to the gap variables here
		for k in range(0, self.__N - 1):
			for j in range(1, len(self.__x[k]) - 1):
				g_pos = k * j * bits
				for i in range(0, bits - 1):
					h = offset + g_pos + i
					self.__bvm[h, h] += self.__l1

	def __addE2Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E1"""
		i = 1
		# TODO Implement E2 coefficients

	def createBVMatrix(self):
		"""Creates the symmetric binary variable matrix of coefficients from the energy function"""
		size = self.calc_maxBV()
		self.__bvm = numpy.zeros((size, size))
		self.__addE0Coefficients()
		self.__addE1Coefficients()
		self.__addE2Coefficients()
		return self.__bvm + self.__bvm.T - numpy.diag(self.__bvm.diagonal())

	def writeQUBO(self, outfilepath, infilepath):
		"""Outputs QUBO format representation"""
		o = open(outfilepath, 'w')

		o.write("c\n")
		o.write("c QUBO format representation of " + infilepath + "\n")
		o.write("c\n")

		# Assume unconstrained target "0"
		o.write("p qubo 0 " + str(self.get_NbBV()) + " " + str(self.get_NbBV()) + " X\n")

		# Output diagonals
		o.write("c\n")
		o.write("c diagonals\n")
		for i in range(0, self.get_NbBV() - 1):
			v = self.__bvm[i, i]
			if v != 0:
				o.write(str(i) + " " + str(i) + " " + str(v) + "\n")

		# Output elements
		o.write("c\n")
		o.write("c elements\n")
		for i in range(0, self.get_NbBV() - 2):
			for j in range(i + 1, self.get_NbBV() - 2):
				v = self.__bvm[i, j]
				if v != 0:
					o.write(str(i) + " " + str(j) + " " + str(v) + "\n")


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
		for i in range(0, self.__bvc.m()):
			val += (self.__bv[offset + i] * 2 ** i)
		return val

	def __G(self, k, j):
		val = 0
		offset = self.__bvc.get_gVarOffset() + (k * self.__bvc.K() + j) * self.__bvc.p()
		for i in range(0, self.__bvc.p()):
			val += (self.__bv[offset + i] * 2 ** i)
		return val

	def __calcE0(self):
		"""Calculates E0 from a defined set of binary variables"""
		x = self.__bvc.x()
		sum_k = 0
		for k in range(0, len(x) - 1):
			sum_j = 0
			for j in range(0, len(x[k])):
				term_1 = self.__x(k, j + 1) - self.__x(k, j) - 1 - self.__G(k, j + 1)
				term_2 = self.__x(k, 0) - self.__G(k, 0)
				sum_j += (term_1 ** 2) + (term_2 ** 2)

			sum_k += sum_j

		return self.__bvc.l0() * sum_k

	def __calcE1(self):
		"""Calculates E1 from a defined set of binary variables"""
		x = self.__bvc.x()
		sum_k = 0.0
		for k in range(0, self.__bvc.N() - 1):
			sum_j = 0.0
			for j in range(1, len(x[k]) - 1):
				sum_j += self.__G(k, j) ** 2
			sum_k += sum_j

		return self.__bvc.l1() * sum_k

	def __calcE2(self):
		"""Calculates E2 from a defined set of binary variables"""
		x = self.__bvc.x()
		N = self.__bvc.N()
		sum_k = 0.0
		for k in range(0, N - 2):
			sum_q = 0.0
			for q in range(k + 1, N - 1):
				sum_i = 0.0
				for i in range(0, len(x[k]) - 1):
					sum_j = 0.0
					for j in range(0, len(x[q]) - 1):
						# sum_j += W[i][j][k][q] * r[i][j][k][q] * math.floor(1 - self.__bvc.l2() * ((self.__x(k,i) - self.__x(q,j)) ** 2))
						# TODO Ignore weighting and reward matrix for now (this obviously won't work in practice!!)
						sum_j += math.floor(1 - self.__bvc.l2() * ((self.__x(k, i) - self.__x(q, j)) ** 2))
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
		for i in range(0, bvlen - 1):
			self.__bv[i] = int(bvstr[i])

		# Value for next iteration
		self.__index += 1

	def calcSolution(self):
		"""Calculates optimal solution by trying all possible permutations of binary variables"""
		min_solution = copy.deepcopy(self.__bv)
		min_val = sys.float_info.max
		nbPermutations = 2 ** len(self.__bv)
		for i in range(1, nbPermutations):
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

	print("Loading input into memory...", end="")
	handle = open(args.input, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	print(" done")
	print()
	print("Input:")

	# Collect variables
	bvc = BVC(args.P, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting)
	for r in records:
		bvc.add_record(r)

		# Input will only be small so just output every sequence found
		print(">" + r.id)
		print(r.seq)

	print()

	# Print current state of variables
	print(bvc)

	# Create matrix
	print()
	print("Creating matrix of coefficients of binary variables ...", end="")
	sys.stdout.flush()
	m = bvc.createBVMatrix()
	# print (m)
	print(" done")

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

			print("Min value is: " + str(val))
			print("Binary variables are:\n" + bv)


main()
