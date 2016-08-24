#!/usr/bin/env python3

"""bvc.py: Manages binary and integer variables for the MSA 2 QUBO problem"""

import math
import numpy as np

try:
    import gurobipy
except ImportError:
	pass


__author__ = "Dan Mapleson, Luis Yanes, Katie Barr, Sophie Kirkwood and Tim Stitt"
__copyright__ = "Copyright 2016, Quantum MSA"
__credits__ = ["Dan Mapleson", "Luis Yanes", "Katie Barr",
                    "Sophie Kirkwood", "Tim Stitt"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"
__status__ = "Prototype"


class BVC:
	"""This class represents the binary variables required for a given number of sequences.  Little endian is used where
	binary variables are representing integer values.
	To use the class initialise, and optionally override default parameters representing the max gap size and weightings,
	then add BioPython sequence records that will represent the MSA.  Finally call createBVMatrix to create the symettric
	matrix containing binary variable coefficients."""


	def __init__(self, P=1, d=2.0, l0=0.8, l1=1.0, l2=10.0, settings_file=None):
		self.data = []

		if settings_file == None:
			self.__P = P
			self.__d = d
			self.__l0 = l0
			self.__l1 = l1
			self.__l2 = l2
			self.__N = 0
			self.__Lmax = 0
			self.__K = 0

			self.__x = []
		else:
			handle = open(settings_file, "r");

			parts = handle.readline().split("\t")
			self.__P = int(parts[1])
			parts = handle.readline().split("\t")
			self.__d = float(parts[1])
			parts = handle.readline().split("\t")
			self.__l0 = float(parts[1])
			parts = handle.readline().split("\t")
			self.__l1 = float(parts[1])
			parts = handle.readline().split("\t")
			self.__l2 = float(parts[1])
			parts = handle.readline().split("\t")
			self.__N = int(parts[1])
			parts = handle.readline().split("\t")
			self.__Lmax = int(parts[1])
			parts = handle.readline().split("\t")
			self.__K = int(parts[1])
			self.__x = []

			for line in handle:
				parts = line.split("\t")
				self.__x.append([0] * int(parts[1]))

			handle.close()


		self.__W = np.zeros((0, 0))

		self.__bvm = np.zeros((0, 0))
		self.__bvmt = np.zeros((0, 0))
		self.__qim = np.zeros((0, 0))
		self.__lil = []

		self.__bvs = []
		self.energy = 0
		self.ienergy = 0
		self.active = []
		self.nb_active = 0

	def __str__(self):
		return 	"N=" + str(self.__N) + "\tNumber of sequences\n" \
				"Lmax=" + str(self.__Lmax) + "\tLength of longest sequence (rounded up to nearest ^2)\n" \
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
				"Max Binary Variables required for rewards: " + str(self.get_NbRewardVars()) + "\n" \
				"Max Binary Variables required for y and z (each): " + str(self.get_NbYVars()) + "\n\n"

	def add_record(self, record):
		"""Adds a biopython record and updates the state of this object accordingly"""
		l = len(record.seq)
		self.__Lmax = max(self.__Lmax, l)
		# ensure Lmax is power of 2
		self.__Lmax = 1<<(self.__Lmax - 1).bit_length()

		self.__x.append([0] * l)  # Add a list the size of the this sequence to the 'x' list
		self.__K += l
		self.__N += 1

	def save_settings(self, filename):
		handle = open(filename, "w");

		print("P\t" + str(self.__P), file=handle)
		print("d\t" + str(self.__d), file=handle)
		print("l0\t" + str(self.__l0), file=handle)
		print("l1\t" + str(self.__l1), file=handle)
		print("l2\t" + str(self.__l2), file=handle)
		print("N\t" + str(self.__N), file=handle)
		print("Lmax\t" + str(self.__Lmax), file=handle)
		print("K\t" + str(self.__K), file=handle)
		for v in self.__x:
			print("x\t" + str(len(v)), file=handle)

		handle.close()

	def load_bvs(self, filename, active):

		with open(filename) as f:
			f.readline()
			s = f.readline().strip()
			#index = 0
			for c in s.strip():
				#while index < len(active) and active[index] == False:
				#	index += 1
				#	self.__bvs.append(0)
				self.__bvs.append(int(c))
				#index += 1

	def get_energy_from_file(self, filename):
		with open(filename) as f:
			f.readline()	# Number of bits
			f.readline()	# Solution
			s = f.readline().strip()	# Energy
			parts = s.split()
			return float(parts[0])


	def make_msa(self):
		msa = []
		b_k = 0
		for k in range(self.__N):
			L_k = len(self.__x[k])
			sa = []
			for j in range(L_k):
				b_kj = b_k + (j * self.m())
				#pos = j # Ensures minimal offset to ensure characters don't stack up on each other (handles the -1 in E0).
				pos = 0
				for a in range(self.m()):
					pos += 2 ** a if self.__bvs[b_kj + a] == 1 else 0
				sa.append(abs(pos))
			msa.append(sa)
			b_k += L_k * self.m()
		return msa

	def make_gap_matrix(self):
		gm = []
		g_k = self.get_gVarOffset()
		for k in range(self.__N):
			klen = len(self.__x[k])
			gr = []
			for j in range(klen):
				g_kj = g_k + (j * self.p())
				pos = 0
				for g_a in range(self.p()):
					pos += 2 ** g_a if self.__bvs[g_kj + g_a] == 1 else 0
				gr.append(pos)
			gm.append(gr)
			g_k += klen * self.p()
		return gm

	def W(self):
		return self.__W

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
		min_m = 2 * self.__Lmax + 1
		return 1<<(min_m - 1).bit_length()

	def Lmax(self):
		return self.__Lmax

	def m(self):
		"""Returns the number of binary variables required to represent each row in solution space"""
		return math.ceil(np.log2(self.M()))

	def p(self):
		"""Returns the number of binary variables required to represent each gap"""
		return math.ceil(np.log2(self.__P + 1))

	def x(self):
		return self.__x

	def l0(self):
		""" Penalty for placing the characters out of order."""
		return self.__l0

	def l1(self):
		""" Penalty for adding gaps."""
		return self.__l1

	def l2(self):
		""" Affinity of characters across sequences, lower values relaxes the aligner to spread-out positions. """
		return self.__l2

	def im(self):
		return self.__qim

	def getPosSolution(self):
		return self.__bvs[0:self.get_NbPositioningVars()]

	def getGapSolution(self):
		return self.__bvs[self.get_gVarOffset():self.get_gVarOffset()+self.get_NbGapVars()]

	def calc_solutionSpaceSize(self):
		"""Size of the solution space"""
		return self.__N * self.M()

	def calc_minBV(self):
		"""Calculate lower bound on number of binary variables required"""
		return self.K() * (self.m() + self.p())

	def calc_maxBV(self):
		"""Calculate upper bound on number of binary variables required"""
		return self.calc_minBV() + math.ceil(self.__N * (self.__N - 1) * (self.__Lmax ** 2) * (1 + 2 * self.m()) / 2)

	def get_NbIV(self):
		return self.get_NbPositioningVars(intmode=True) + self.get_NbGapVars(intmode=True)

	def get_NbBV(self):
		"""Return the actual number of binary variables required for this problem"""
		#return self.get_NbPositioningVars() + self.get_NbGapVars() + self.get_NbRewardVars() + self.get_NbYVars() + self.get_NbZVars()
		return self.get_NbPositioningVars() + self.get_NbGapVars()

	def get_gVarOffset(self, intmode=False):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a gap"""
		return self.get_NbPositioningVars(intmode)

	def get_rVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_gVarOffset() + self.get_NbGapVars()

	def get_yVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_rVarOffset() + self.get_NbRewardVars()

	def get_zVarOffset(self):
		"""Gets the offset (in terms of number of binary variables) for the first binary variable representing a reward"""
		return self.get_yVarOffset() + self.get_NbYVars()

	def get_NbPositioningVars(self, intmode=False):
		return self.K() * self.m() if not intmode else self.K()

	def get_NbGapVars(self, intmode=False):
		return self.K() * self.p() if not intmode else self.K()

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

	def lenK(self, k):
		return len(self.__x[k])

	def qim(self, x, y):
		return self.__qim[x, y]

	def lil(self, x):
		return self.__lil[x]

	def printqim(self):
		for k in range(self.N()):
			L_k = len(self.__x[k])
			for j in range(L_k):
				print("x_", k, ",", j, sep="", end="")
				print(" ", end="")
		for k in range(self.N()):
			L_k = len(self.__x[k])
			for j in range(L_k):
				print("G_", k, ",", j, sep="", end="")
				print(" ", end="")
		print()
		print(self.__qim)

	def printlil(self):
		for k in range(self.N()):
			L_k = len(self.__x[k])
			for j in range(L_k):
				print("x_", k, ",", j, sep="", end="")
				print(" ", end="")
		for k in range(self.N()):
			L_k = len(self.__x[k])
			for j in range(L_k):
				print("G_", k, ",", j, sep="", end="")
				print(" ", end="")
		print()
		print(self.__lil)

	def printIntegerCoefficients(self):
		np.set_printoptions(threshold=np.inf, linewidth=np.inf, suppress=True, precision=2)
		print("Quadratic coefficients:")
		self.printqim()
		print()
		print("Linear coefficients:")
		self.printlil()
		print()
		print("Constant:")
		print(self.ienergy)
		print()


	def __addE0Coefficients(self, intmode=False):
		"""Updates the binary variable matrix with coefficients from E0"""

		if intmode:
			x_k = 0
			G_k = self.get_gVarOffset(True)
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(L_k - 1):
					x_kj = x_k + j
					x_kj1 = x_kj + 1
					G_kj1 = G_k + j + 1

					# Nodes
					self.__qim[x_kj, x_kj] += 1		# Squared
					self.__qim[x_kj1, x_kj1] += 1	# Squared
					self.__qim[G_kj1, G_kj1] += 1	# Squared
					self.__qim[x_k, x_k] += 1		# Squared
					self.__qim[G_k, G_k] += 1		# Squared

					# Couplers
					self.__qim[x_kj, x_kj1] += -2
					self.__qim[x_kj, G_kj1] += 2
					self.__qim[x_kj1, G_kj1] += -2
					self.__qim[x_k, G_k] += -2

					# Linear parts
					self.__lil[x_kj] += 2
					self.__lil[x_kj1] -= 2
					self.__lil[G_kj1] += 2

					# Left over energy
					self.ienergy += 1

				x_k += L_k
				G_k += L_k

			nbvars = self.get_NbPositioningVars(intmode=True) + self.get_NbGapVars(intmode=True)
			for k in range(nbvars):
				for j in range(k, nbvars):
					self.__qim[k, j] *= self.__l0
				self.__lil[k] *= self.__l0
			self.ienergy *= self.__l0
		else:

			x_k = 0
			G_k = self.get_gVarOffset()
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(L_k - 1):
					x_kj = x_k + (j * self.m())
					x_kj1 = x_k + ((j + 1) * self.m())
					G_kj1 = G_k + ((j + 1) * self.p())
					self.energy += 1

					# x nodes.  Include quadratic and linear terms
					for x_a1 in range(self.m()):
						for x_a2 in range(x_a1, self.m()):
							x_kja1 = x_kj + x_a1
							x_kj1a1 = x_kj1 + x_a1
							x_kja2 = x_kj + x_a2
							x_kj1a2 = x_kj1 + x_a2
							x_k1a1 = x_k + x_a1
							x_k1a2 = x_k + x_a2

							quad_scale = (2 ** x_a1) ** 2 if x_a1 == x_a2 else (2 ** (x_a1 + x_a2 + 1))

							# Quadratic parts
							self.__bvm[x_kja1, x_kja2] += quad_scale
							self.__bvm[x_kj1a1, x_kj1a2] += quad_scale
							self.__bvm[x_k1a1, x_k1a2] += quad_scale

							# Linear parts
							if x_a1 == x_a2:
								self.__bvm[x_kja1, x_kja2] += 2 * (2 ** x_a1)
								self.__bvm[x_kj1a1, x_kj1a2] -= 2 * (2 ** x_a1)

					# g node
					for G_a1 in range(self.p()):
						for G_a2 in range(G_a1, self.p()):
							g_kj1a1 = G_kj1 + G_a1
							g_kj1a2 = G_kj1 + G_a2
							g_k1a1 = G_k + G_a1
							g_k1a2 = G_k + G_a2

							quad_scale = (2 ** G_a1) ** 2 if G_a1 == G_a2 else (2 ** (G_a1 + G_a2 + 1))

							# Quadratic parts
							self.__bvm[g_kj1a1, g_kj1a2] += quad_scale
							self.__bvm[g_k1a1, g_k1a2] += quad_scale

							# Linear parts
							if G_a1 == G_a2:
								self.__bvm[g_kj1a1, g_kj1a2] += 2 * (2 ** G_a1)

					# xx - coupled
					for x_a1 in range(self.m()):
						for x_a2 in range(self.m()):
							quad_scale = (2 ** x_a1) ** 2 if x_a1 == x_a2 else (2 ** (x_a1 + x_a2))
							self.__bvm[x_kj + x_a1, x_kj1 + x_a2] -= 2 * quad_scale


					# xG - coupled
					for x_a in range(self.m()):
						for G_a in range(x_a, self.p()):
							quad_scale = (2 ** x_a) ** 2 if x_a == G_a else (2 ** (G_a + x_a))
							self.__bvm[x_kj + x_a, G_kj1 + G_a] += 2 * quad_scale
							self.__bvm[x_kj1 + x_a, G_kj1 + G_a] -= 2 * quad_scale
							self.__bvm[x_k + x_a, G_k + G_a] -= 2 * quad_scale

				x_k += L_k * self.m()
				G_k += L_k * self.p()

			nbvars = self.get_NbPositioningVars() + self.get_NbGapVars()
			for k in range(nbvars):
				for j in range(k, nbvars):
					self.__bvm[k, j] *= self.__l0
			self.energy *= self.__l0




	def __addE1Coefficients(self, intmode=False):
		"""Updates the binary variable matrix with coefficients from E1"""

		if intmode:
			g_k = self.get_gVarOffset(True)
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(1, L_k):
					g_kj = g_k + j
					self.__qim[g_kj, g_kj] += 1
				g_k += L_k

			nbvars = self.get_NbPositioningVars(intmode=True) + self.get_NbGapVars(intmode=True)
			for k in range(self.get_NbPositioningVars(intmode=True)+1,nbvars):
				for j in range(k, nbvars):
					self.__qim[k, j] *= self.__l1

		else:
			g_k = self.get_gVarOffset()
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(1, L_k):
					g_kj = g_k + (j * self.p())
					for a1 in range(self.p()):
						for a2 in range(a1, self.p()):
							g_kja1 = g_kj + a1
							g_kja2 = g_kj + a2
							quad_scale = (2 ** a1) ** 2 if a1 == a2 else (2 ** (a1 + a2 + 1))
							self.__bvm[g_kja1, g_kja2] += self.__l1*quad_scale
				g_k += L_k * self.p()

	def __addE2Coefficients(self, intmode=False):
		"""Updates the binary variable matrix with coefficients from E2"""
		e2m = np.zeros
		if intmode:
			pass
		else:
			pass

	def __calcNbNodes(self):
		count = 0
		for i in range(self.get_NbBV()):
			if self.__bvm[i, i] != 0.0:
				count += 1
		return count

	def __calcNbCouplers(self):
		count = 0
		for i in range(self.get_NbBV() - 1):
			for j in range(i + 1, self.get_NbBV()):
				if self.__bvm[i, j] != 0.0:
					count += 1
		return count

	def calcActiveBVs(self):
		self.active = [False] * self.get_NbBV()
		for i in range(self.get_NbBV()):
			for j in range(i, self.get_NbBV()):
				if self.__bvm[i, j] != 0.0:
					self.active[i] = True
					break
		self.nb_active = sum(self.active)
		return self.nb_active

	def createW(self):

		N=self.__N
		self.__W=np.zeros(shape=(self.__Lmax,self.__Lmax,N,N))

		for k in range(N-1):
			for q in range(k+1,N):
				for i in range(len(self.__x[k])):
					for j in range(len(self.__x[q])):
						self.__W[i,j,k,q] = 1. if self.__x[k][i] == self.__x[q][j] else 0.

	def w(self, i, j, k, q):
		return self.__W[i,j,k,q]

	def sophiesMethod(self):
		bvec = np.zeros((1, len(self.__bvm))) #1D array
		np.put(bvec, [x for x in range(0, len(self.__bvm))], [0, 0, 0, 0, 0,
															  1, 0, 0, 0, 0,
															  0, 1, 0, 0, 0,
															  1, 1, 0, 0, 0,
															  0, 0, 1, 0, 0,
															  0, 0, 0, 0, 0])
		print(bvec)
		#np.put(bvec, ['''index values'''], ['''vector values'''])
		print("Optimal solution")
		print(np.dot(bvec, np.dot(self.__bvm, bvec.transpose())))

	def createBVMatrix(self, intmode=False):
		"""Creates the symmetric binary variable matrix of coefficients from the energy function"""
		if not intmode:
			size = self.get_NbBV()
			self.__bvm = np.zeros((size, size))
			self.__bvmt = np.zeros((size, size))
			self.energy = 0

			np.set_printoptions(threshold=np.inf, linewidth=np.inf, suppress=True, precision=2)

			self.__addE0Coefficients()
			#print("\n\nBVM - After E0 applied\n", self.__bvm)
			self.__addE1Coefficients()
			#print("\n\nBVM - After E1 applied\n", self.__bvm)
			#self.__addE2Coefficients()

			#self.sophiesMethod()

		else:
			isize = self.get_NbIV()
			self.__qim = np.zeros((isize, isize))
			self.__lil = [0] * isize
			self.ienergy = 0
			self.__addE0Coefficients(intmode=True)
			self.__addE1Coefficients(intmode=True)

		#return self.__bvm + self.__bvm.T - numpy.diag(self.__bvm.diagonal())

	def writeQUBO(self, outfilepath, infilepath):
		"""Outputs QUBO format representation"""
		o = open(outfilepath, 'w')

		print("c ---------\nc", file=o)
		print("c QUBO format representation of ", infilepath, file=o)
		print("c\np qubo 0 ", self.get_NbBV(), self.__calcNbNodes(), self.__calcNbCouplers(), file=o)
		print("c\nc nodes\nc", file=o)
		for i in range(self.get_NbBV()):
			v = self.__bvm[i, i]
			if v != 0:
				print(i, i, v, file=o)

		# Output elements
		print("c\nc couplers\nc", file=o)
		for i in range(self.get_NbBV() - 1):
			for j in range(i + 1, self.get_NbBV()):
				v = self.__bvm[i, j]
				if v != 0:
					print(i, j, v, file=o)

		o.close()

