#!/usr/bin/env python3

import os
import argparse
import math
import numpy
import copy
import sys
from Bio import SeqIO


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


		self.__W = numpy.zeros((0, 0))

		self.__bvm = numpy.zeros((0, 0))
		self.__im = numpy.zeros((0, 0))

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

	def im(self):
		return self.__im

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

	def get_NbBV(self):
		"""Return the actual number of binary variables required for this problem"""
		return self.get_NbPositioningVars() + self.get_NbGapVars() + self.get_NbRewardVars() + self.get_NbYVars() + self.get_NbZVars()

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

	def get_NbGapVars(self):
		return self.K() * self.p()

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
					self.__im[x_kj, x_kj] += self.__l0 ** 2 + 2 * self.__l0
					self.__im[x_kj1, x_kj1] += self.__l0 ** 2 - 2 * self.__l0
					self.__im[G_kj1, G_kj1] += self.__l0 ** 2 + self.__l0
					self.__im[x_kj, x_kj1] -= 2 * self.__l0 * self.__l0
					self.__im[x_kj, G_kj1] += 2 * self.__l0 * self.__l0
					self.__im[x_kj1, G_kj1] -= 2 * self.__l0 * self.__l0
					self.ienergy += 1.0 * self.__l0
					self.__im[x_k, x_k] += self.__l0 ** 2
					self.__im[G_k, G_k] += self.__l0 ** 2
					#self.__bvm[b_ka, g_ka] -= 2 * v_b_a * v_g_a
				x_k += L_k
				G_k += L_k
		else:
			x_k = 0
			G_k = self.get_gVarOffset()
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(L_k - 1):
					x_kj = x_k + (j * self.m())
					x_kj1 = x_k + ((j + 1) * self.m())
					G_kj1 = G_k + ((j + 1) * self.p())
					self.energy += 1.0 * self.__l0

					# x node
					for x_a in range(self.m()):
						x_k1a = x_k + x_a
						x_kja = x_kj + x_a
						x_kj1a = x_kj1 + x_a
						x_scale = 2 ** x_a
						self.__bvm[x_kja, x_kja] += (x_scale ** 2 + 2 * x_scale) * self.__l0
						self.__bvm[x_kj1a, x_kj1a] += (x_scale ** 2 - 2 * x_scale) * self.__l0
						self.__bvm[x_k1a, x_k1a] += (x_scale ** 2) * self.__l0

					# g node
					for G_a in range(self.p()):
						g_k1a = G_k + G_a
						g_kj1a = G_kj1 + G_a
						g_scale = 2 ** G_a
						self.__bvm[g_kj1a, g_kj1a] += (g_scale ** 2 + 2 * g_scale) * self.__l0
						self.__bvm[g_k1a, g_k1a] += g_scale * self.__l0

					# xx - coupled
					for x_a1 in range(self.m()):
						x_scale1 = 2 ** x_a1
						for x_a2 in range(self.m()):
							x_scale2 = 2 ** x_a2
							self.__bvm[x_kj + x_a1, x_kj + x_a2] -=  2 * x_scale1 * x_scale2 * self.__l0

					# xy - coupled
					for x_a in range(self.m()):
						x_scale = self.__l0 * 2 ** x_a
						for G_a in range(self.p()):
							g_scale = 2 ** G_a
							self.__bvm[x_kj + x_a, G_kj1 + G_a] += 2 * x_scale * g_scale * self.__l0
							self.__bvm[x_kj1 + x_a, G_kj1 + G_a] -= 2 * x_scale * g_scale * self.__l0
							self.__bvm[x_k + x_a, G_k + G_a] -= 2 * x_scale * g_scale * self.__l0

				x_k += L_k * self.m()
				G_k += L_k * self.p()


	def __addE1Coefficients(self, intmode=False):
		"""Updates the binary variable matrix with coefficients from E1"""

		if intmode:
			g_k = self.get_gVarOffset(True)
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(1, L_k):
					g_kj = g_k + j
					self.__im[g_kj, g_kj] += self.__l1
				g_k += L_k
		else:
			g_k = self.get_gVarOffset()
			for k in range(self.__N):
				L_k = len(self.__x[k])
				for j in range(1, L_k):
					g_kj = g_k + (j * self.p())
					for a in range(self.p()):
						g_kja = g_kj + a
						self.__bvm[g_kja, g_kja] += (2 ** a) * self.__l1
				g_k += L_k * self.p()

	def __addE2Coefficients(self):
		"""Updates the binary variable matrix with coefficients from E2"""
		for k in range(self.__N - 1):
			L_k = len(self.__x[k])
			k_offset = k * (self.__N - 1) * self.__Lmax ** 2
			r_k = self.get_rVarOffset() + k_offset
			y_k = self.get_yVarOffset() + (k_offset * self.m())
			z_k = self.get_zVarOffset() + (k_offset * self.m())
			for q in range(k + 1, self.__N):
				qlen = len(self.__x[q])
				b_k = k * L_k * self.m()
				b_q = q * qlen * self.m()
				e = (q - k - 1) * self.__Lmax ** 2
				r_kq = r_k + e
				y_kq = y_k + (e * self.m())
				z_kq = z_k + (e * self.m())
				for i in range(L_k):
					r_ikq = r_kq + (i * self.__Lmax)
					y_ikq = y_kq + (i * self.__Lmax * self.m())
					z_ikq = z_kq + (i * self.__Lmax * self.m())
					for j in range(qlen):
						b_ki = b_k + (i * self.m())
						b_qj = b_q + (j * self.m())
						w_ijkq = self.__W[i,j,k,q]
						wl2 = w_ijkq * self.__l2
						wl2d = wl2 * self.__d
						r_ijkq = r_ikq + j
						y_ijkq = y_ikq + (j * self.m())
						z_ijkq = z_ikq + (j * self.m())
						for b_a in range(self.m()):
							b_kia = b_ki + b_a
							b_qja = b_qj + b_a
							y_ijkqa = y_ijkq + b_a
							z_ijkqa = z_ijkq + b_a

							self.__bvm[r_ijkq, r_ijkq] += w_ijkq

							self.__bvm[y_ijkqa, b_kia] -= wl2
							self.__bvm[y_ijkqa, y_ijkqa] -= 3 * wl2d
							self.__bvm[r_ijkq, b_kia] -= wl2d
							self.__bvm[y_ijkqa, r_ijkq] += 2 * wl2d
							self.__bvm[y_ijkqa, b_kia] += 2 * wl2d

							self.__bvm[z_ijkqa, b_qja] += wl2
							self.__bvm[z_ijkqa, z_ijkqa] -= 3 * wl2d
							self.__bvm[r_ijkq, b_qja] -= wl2d
							self.__bvm[z_ijkqa, r_ijkq] += 2 * wl2d
							self.__bvm[z_ijkqa, b_qja] += 2 * wl2d

							self.__bvm[y_ijkqa, b_qja] -= 2 * wl2
							self.__bvm[y_ijkqa, y_ijkqa] -= 3 * wl2d
							self.__bvm[r_ijkq, b_kia] -= wl2d
							self.__bvm[y_ijkqa, r_ijkq] += 2 * wl2d
							self.__bvm[y_ijkqa, b_kia] += 2 * wl2d


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
		self.__W=numpy.empty(shape=(self.__Lmax,self.__Lmax,N,N))

		for k in range(N-1):
			for q in range(k+1,N):
				for i in range(len(self.__x[k])):
					for j in range(len(self.__x[q])):
						self.__W[i,j,k,q] = 1 if self.__x[k][i] == self.__x[q][j] else 0

	def w(self, i, j, k, q):
		return self.__W[i,j,k,q]

	def createBVMatrix(self):
		"""Creates the symmetric binary variable matrix of coefficients from the energy function"""
		size = self.get_NbBV()
		self.__bvm = numpy.zeros((size, size))
		self.__addE0Coefficients()
		self.__addE1Coefficients()
		#self.__addE2Coefficients()
		if True:
			isize = self.K() * 2
			self.__im = numpy.zeros((isize, isize))
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
