#!/usr/bin/env python3

import os
import argparse
import math
import numpy
import copy
import sys
from Bio import SeqIO
from bvc import BVC


def main():
	parser = argparse.ArgumentParser("Convert QUBO output (for now a list of line separated 0's and 1's) into an MSA")
	parser.add_argument("settings", help="The file containing settings used to generate the QUBO problem")
	parser.add_argument("bv", help="The file containing solved QUBO binary variables")
	parser.add_argument("input", help="The original input file in Fasta format")
	parser.add_argument("-o", "--output", help="The output file, containing the solved MSA")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	bvc = BVC(settings_file=args.settings)
	print("Loaded QUBO settings")
	print()
	bvc.load_bvs(args.bv)
	print("Loaded binary variables")

	print("Loading input into memory...", end="")
	handle = open(args.input, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	print(" done")
	print()

	print("Settings:")
	print(bvc)
	print()

	msa = bvc.make_msa()
	print("Made MSA")

	print()
	print("Position matrix:")
	for sa in msa:
		for base in sa:
			print(base, end="")
		print()

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

main()
