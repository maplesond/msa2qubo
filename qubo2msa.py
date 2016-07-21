#!/usr/bin/env python3

import os
import argparse
import math
import numpy
import copy
import sys
from Bio import SeqIO
from bvc import BVC

class Qubo2Msa:

	def __init__(self, settings, solution, input, output, active, verbose):
		self.data = []
		self.settings = settings
		self.solution = solution
		self.input = input
		self.output = output
		self.active = active
		self.verbose = verbose

	def run(self):

		bvc = BVC(settings_file=self.settings)
		print("Loaded QUBO settings")
		print()
		bvc.load_bvs(self.solution, self.active)
		print("Loaded solution to QUBO problem with", len(self.active), "binary variables (", sum(self.active), "of which are active. )")
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
