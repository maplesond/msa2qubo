#!/usr/bin/env python3

"""qubo2msa.py: Converts output from a D-Wave device or D-Wave simulator to a MSA."""

import argparse
import numpy as np
from Bio import SeqIO
from bvc import BVC

__author__ = "Dan Mapleson, Luis Yanes, Katie Barr, Sophie Kirkwood and Tim Stitt"
__copyright__ = "Copyright 2016, Quantum MSA"
__credits__ = ["Dan Mapleson", "Luis Yanes", "Katie Barr",
                    "Sophie Kirkwood", "Tim Stitt"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"
__status__ = "Prototype"

class Qubo2Msa:

	def __init__(self, settings, solution, input, output, active, target_energy, verbose):
		self.data = []
		self.settings = settings
		self.solution = solution
		self.input = input
		self.output = output
		self.active = active
		self.verbose = verbose
		self.target_energy = target_energy
		self.otherbvm = None

	def bvm(self, bvc):
		self.otherbvm = bvc

	def run(self):

		bvc = BVC(settings_file=self.settings)
		print("Loaded QUBO settings")
		print()
		bvc.load_bvs(self.solution, self.active)
		print("Loaded solution to QUBO problem with", len(self.active), "binary variables (", sum(self.active), "of which are active. )")
		print()
		energy = bvc.get_energy_from_file(self.solution)
		print("Energy - Target:", self.target_energy)
		print("Energy - Actual:", energy)
		print("Tolerance: 0.5")
		diff = abs(energy - self.target_energy)
		if diff < 0.5:
			print("Difference between actual and target energy is within tolerance:", diff)
		else:
			print("*****************************************************************************")
			print("WARNING: Difference between actual and target energy exceeds tolerance:", diff)
			print("WARNING: It is likely that the solution provided by the solver is not optimal")
			print("*****************************************************************************")

		print()
		print("Loading input sequences into memory...", end="")
		handle = open(self.input, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		print(" done")
		print()

		if self.verbose:
			print("Settings:")
			print(bvc)
			print()

		msa = bvc.make_msa()
		print("Made MSA")

		if self.verbose:
			print()
			print()
			print("Solution variables:")
			print(bvc.getSolutionVars())

			print("Solution shape:")
			print(bvc.getSolutionShape())
			print()
			print()

			if self.otherbvm:
				x=np.reshape(np.asarray(bvc.getSolutionVars()), newshape=(1,len(bvc.getSolutionVars())))
				self.otherbvm.sophiesMethod(x)

		print("Position variables:")
		print(bvc.getPosSolution())
		print()
		print("Position matrix:")
		for sa in msa:
			print(sa)

		print()
		print("Gap variables:")
		print(bvc.getGapSolution())

		gm = bvc.make_gap_matrix()
		print()
		print("Gap matrix:")
		for g in gm:
			print(g)

		# Shouldn't need this when run on the real thing but for now on random data M is not necesarily sufficient to
		# hold the potential position values
		width = 2 ** bvc.m()

		print()
		print("MSA:")
		ss = [[" " for x in range(width)] for x in range(bvc.N())]

		for i in range(bvc.N()):
			sa = msa[i]
			rec = records[i]
			if not len(sa) == len(rec):
				print("ERROR")
				exit(1)

			for j in range(len(sa)):
				pos = sa[j]
				base = rec.seq[j]
				ss[i][pos] = base

		for i in range(bvc.N()):
			for j in range(width):
				print(ss[i][j], end="")
			print()



def main():
	parser = argparse.ArgumentParser("Convert QUBO output (for now a list of line separated 0's and 1's) into an MSA")
	parser.add_argument("settings", help="The file containing settings used to generate the QUBO problem")
	parser.add_argument("bv", help="The file containing solved QUBO binary variables")
	parser.add_argument("input", help="The original input file in Fasta format")
	parser.add_argument("-o", "--output", help="The output file, containing the solved MSA")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	q2m = Qubo2Msa(settings = args.settings, solution=args.bv, input=args.input, output=args.output, verbose=args.verbose)
	q2m.run()

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
