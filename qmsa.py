#!/usr/bin/env python3

import os
import argparse
import subprocess
from msa2qubo import Msa2Qubo
from qubo2msa import Qubo2Msa

def main():
	parser = argparse.ArgumentParser("Perform MSA using D-Wave Simulator")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output_dir", required=True,
						help="The output directory")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)

	print("Converting problem into QUBO form")
	print("---------------------------------")
	print()

	m2q = Msa2Qubo(input=args.input, output=args.output_dir + "/qmsa.qubo", P=args.P, verbose=args.verbose, delta=2.0, l0=0.8, l1=0.1, l2=10.0)
	m2q.run()

	print()
	print("Solving QUBO problem")
	print("--------------------")
	print()

	# Assume qbsolv in on the PATH
	subprocess.call(['qbsolv', '-i', args.output_dir + "/qmsa.qubo", "-o", args.output_dir + "/qmsa.solution"])

	print()
	print("Interpretting from solution")
	print("---------------------------")
	print()

	m2q = Qubo2Msa(settings=args.output_dir + "/qmsa.qubo.settings", solution=args.output_dir + "/qmsa.solution", input=args.input, output=args.output_dir + "/qmsa.msa", verbose=args.verbose)
	m2q.run()



main()
