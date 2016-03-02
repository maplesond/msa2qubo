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

	def __init__(self, P=1, d=2.0, l0=0.8, l1=1.0, l2=10.0):
		self.data = []
		self.__P = P
		self.__d = d
		self.__l0 = l0
		self.__l1 = l1
		self.__l2 = l2

	def __str__(self):
		return 	"N=" + str(self.__N) + "\tNumber of sequences\n" \
				"Lmax=" + str(self.__Lmax) + "\tLength of longest sequence\n" \
				"K=" + str(self.__K) + "\tTotal number of characters\n" \
				"M=" + str(self.M()) + "\tLength of each row in solution space\n" \
				"m=" + str(self.m()) + "\tNumber of binary variables required to represent each position in solution space\n" \
				"P=" + str(self.__P) + "\tMaximum gap size\n" \
				"p=" + str(self.p()) + "\tNumber of binary variables required to represent each gap\n\n" \
				"l0=" + str(self.__l0) + "\tPosition weighting\n" \
				"l1=" + str(self.__l1) + "\tGap weighting\n" \
				"l2=" + str(self.__l2) + "\tReward weighting\n" \
				"d=" + str(self.__d) + "\tScaling factor when substituting products of 3BVs to 2BVs\n\n" \
				"Solution space will contain " + str(self.calc_solutionSpaceSize()) + " cells\n\n" \
				"Bounds on total # Binary Variables required: " + str(self.calc_minBV()) + "-" + str(self.calc_maxBV()) + "\n" \
				"# Binary Variables required for positioning: " + str(self.get_NbPositioningVars()) + "\n" \
				"# Binary Variables required for gaps: " + str(self.get_NbGapVars()) + "\n" \
				"# Binary Variables required for rewards: " + str(self.get_NbRewardVars()) + "\n" \
				"# Binary Variables required for y and z (each): " + str(self.get_NbYVars()) + "\n" \
				"Total # Binary Variables required: " + str(self.get_NbBV()) + "\n" \

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

	def d(self):
		return self.__d

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
		return self.get_NbPositioningVars() + self.get_NbGapVars() + self.get_NbRewardVars() + self.get_NbYVars() + self.get_NbZVars()

	def get_gVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a gap"""
		return self.__K * self.m()

	def get_rVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_gVarOffset() + self.get_NbGapVars()

	def get_yVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_rVarOffset() + self.get_NbRewardVars()

	def get_zVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_yVarOffset() + self.get_NbYVars()

	def get_NbPositioningVars(self):
		return self.__K * self.m()

	def get_NbGapVars(self):
		return self.__K * self.p()

	def get_NbRewardVars(self):
		count = 0
		for k in range(self.__N - 1):
			klen = len(self.__x[k])
			for q in range(k + 1, self.__N):
				qlen = len(self.__x[q])
				for i in range(klen):
					for j in range(qlen):
						count += 1
		return count

	def get_NbYVars(self):
		return self.get_NbRewardVars() * self.m()

	def get_NbZVars(self):
		return self.get_NbRewardVars() * self.m()

	def __addE0Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E0"""
		for k in range(self.__N):
			klen = len(self.__x[k])
			b_k = k * klen * self.m()
			g_k = self.get_gVarOffset() + k * klen * self.p()
			for j in range(klen - 1):
				b_kj = b_k + (j * self.m())
				b_kj1 = b_k + ((j + 1) * self.m())
				g_kj1 = g_k + ((j + 1) * self.p())
				for b_a in range(self.m()):
					b_kja = b_kj + b_a
					b_kj1a = b_kj1 + b_a
					for g_a in range(self.p()):
						g_kj1a = g_kj1 + g_a
						self.__bvm[b_kj1a, b_kj1a] += self.__l0 ** 2
						self.__bvm[b_kj1a, b_kja] -= self.__l0
						self.__bvm[b_kj1a, b_kj1a] -= self.__l0
						self.__bvm[b_kj1a, g_kj1a] -= self.__l0
						self.__bvm[b_kja, b_kj1a] -= self.__l0
						self.__bvm[b_kja, b_kja] += self.__l0 ** 2
						self.__bvm[b_kja, b_kja] -= self.__l0
						self.__bvm[b_kja, g_kj1a] += self.__l0
						self.__bvm[b_kj1a, b_kj1a] -= self.__l0
						self.__bvm[b_kja, b_kja] -= self.__l0
						# +1 here??? What do wath a constant?
						self.__bvm[g_kj1a, g_kj1a] += self.__l0
						self.__bvm[g_kj1a, b_kj1a] -= self.__l0
						self.__bvm[g_kj1a, b_kja] += self.__l0
						self.__bvm[g_kj1a, g_kj1a] += self.__l0
						self.__bvm[g_kj1a, g_kj1a] += self.__l0 ** 2

				# Coefficients for second term (gap prior to first character)
				for b_a in range(self.m()):
					for g_a in range(self.p()):
						b_ka = b_k + b_a
						g_ka = g_k + g_a
						self.__bvm[b_ka, b_ka] += self.__l0 ** 2
						self.__bvm[g_ka, g_ka] += self.__l0 ** 2
						self.__bvm[b_ka, g_ka] -= self.__l0
						self.__bvm[g_ka, b_ka] -= self.__l0

	def __addE1Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E1"""
		for k in range(self.__N):
			klen = len(self.__x[k])
			g_k = self.get_gVarOffset() + k * klen * self.p()
			for j in range(1, klen):
				g_kj = g_k + (j * self.p())
				for a in range(0, self.p() - 1):
					g_kja = g_kj + a
					print(g_kja)
					self.__bvm[g_kja, g_kja] += ((2 ** a) * self.__l1) ** 2

	def __addE2Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E2"""
		for k in range(self.__N - 1):
			for q in range(k + 1, self.__N):
				klen = len(self.__x[k])
				b_k = k * klen * self.m()
				b_q = q * klen * self.m()
				for i in range(klen):
					qlen = len(self.__x[q])
					for j in range(qlen):
						b_ki = b_k + (i * self.m())
						b_qj = b_q + (j * self.m())
						r_ijkq = self.get_rVarOffset() + k * j
						for b_a in range(self.m()):
							b_kia = b_ki + b_a
							b_qja = b_qj + b_a
							#y_ijkqa = self.get_yVarOffset() +
							#z_ijkqa = self.get_zVarOffset()
							#self.__bvm[h, h] += self.__l1

						# sum_j += W[i][j][k][q] * r[i][j][k][q] * math.floor(1 - self.__bvc.l2() * ((self.__x(k,i) - self.__x(q,j)) ** 2))
						# TODO Ignore weighting and reward matrix for now (this obviously won't work in practice!!)
						#sum_j += math.floor(1 - self.__bvc.l2() * ((self.__x(k, i) - self.__x(q, j)) ** 2))
			

	def __calcNbDiagonals(self):
		count = 0
		for i in range(self.get_NbBV()):
			if self.__bvm[i, i] != 0.0:
				count += 1
		return count

	def __calcNbElements(self):
		count = 0
		for i in range(self.get_NbBV() - 1):
			for j in range(i + 1, self.get_NbBV()):
				if self.__bvm[i, j] != 0.0:
					count += 1
		return count

	def createBVMatrix(self):
		"""Creates the symmetric binary variable matrix of coefficients from the energy function"""
		size = self.get_NbBV()
		self.__bvm = numpy.zeros((size, size))
		self.__addE0Coefficients()
		self.__addE1Coefficients()
		self.__addE2Coefficients()
		#return self.__bvm + self.__bvm.T - numpy.diag(self.__bvm.diagonal())

	def writeQUBO(self, outfilepath, infilepath):
		"""Outputs QUBO format representation"""
		o = open(outfilepath, 'w')

		o.write("c\n")
		o.write("c QUBO format representation of " + infilepath + "\n")
		o.write("c\n")

		# Assume unconstrained target "0"
		o.write("p qubo 0 " + str(self.get_NbBV()) + " " + str(self.__calcNbDiagonals()) + " " + str(self.__calcNbElements()) + "\n")

		# Output diagonals
		o.write("c\n")
		o.write("c diagonals\n")
		for i in range(self.get_NbBV()):
			v = self.__bvm[i, i]
			if v != 0:
				o.write(str(i) + " " + str(i) + " " + str(v) + "\n")

		# Output elements
		o.write("c\n")
		o.write("c elements\n")
		for i in range(self.get_NbBV() - 1):
			for j in range(i + 1, self.get_NbBV()):
				v = self.__bvm[i, j]
				if v != 0:
					o.write(str(i) + " " + str(j) + " " + str(v) + "\n")

		o.close()


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
	parser.add_argument("-d", type=float, default=2.0, help="Delta.  The scaling factor to apply when converting products of 3BVs to 2BVs.")
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
	bvc = BVC(P=args.P, d=args.d, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting)
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
	bvc.createBVMatrix()
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
