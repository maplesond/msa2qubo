#!/usr/bin/env python3

import os
import argparse
import subprocess
import time
from msa2qubo import Msa2Qubo
from qubo2msa import Qubo2Msa

def main():
	parser = argparse.ArgumentParser("Perform MSA using D-Wave Simulator")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output_dir", required=True,
						help="The output directory")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA")
	parser.add_argument("-l0", "--position_weighting", type=float, default=1.0,
						help="The weighting to apply to positioning of elements (must be larger than --gap_weighting)")
	parser.add_argument("-l1", "--gap_weighting", type=float, default=0.1,
						help="The weighting to apply to gap penalties")
	parser.add_argument("-l2", "--reward_weighting", type=float, default=10.0,
						help="The weighting to apply to reward matches (must be greater than 1.0)")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	startall = time.time()

	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)


	print("Converting problem into QUBO form")
	print("---------------------------------")
	print()

	start1 = time.time()
	m2q = Msa2Qubo(input=args.input, output=args.output_dir + "/qmsa.qubo", P=args.P, verbose=args.verbose, delta=2.0, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting)
	m2q.run()
	end1 = time.time()
	print("Time taken to create QUBO file (s): ", "{0:.2f}".format(round(end1 - start1,2)))

	print()
	print()
	print("Solving QUBO problem")
	print("--------------------")
	print()

	# Assume qbsolv in on the PATH
	start2 = time.time()
	cmd_args = ['qbsolv', '-i', args.output_dir + "/qmsa.qubo", "-o", args.output_dir + "/qmsa.solution"]
	print("Executing:", " ".join(cmd_args))
	subprocess.call(cmd_args)
	subprocess.call(['cat', args.output_dir + "/qmsa.solution"])
	end2 = time.time()
	print("Time taken to solve QUBO problem (s): ", "{0:.2f}".format(round(end2 - start2,2)))

	print()
	print()
	print("Interpretting from solution")
	print("---------------------------")
	print()

	start3 = time.time()
	m2q = Qubo2Msa(settings=args.output_dir + "/qmsa.qubo.settings", solution=args.output_dir + "/qmsa.solution", input=args.input, output=args.output_dir + "/qmsa.msa", active=m2q.active, verbose=args.verbose)
	m2q.run()
	end3 = time.time()
	print("Time taken to interpret solution (s): ", "{0:.2f}".format(round(end3 - start3,2)))

	endall = time.time()
	print()
	print("Total time taken: ", "{0:.2f}".format(round(endall - startall,2)))

main()