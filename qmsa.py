#!/usr/bin/env python3

"""qmsa.py: Pipeline script for driving MSA via a D-Wave simulator"""

import os
import argparse
import subprocess
import time
from msa2qubo import Msa2Qubo
from qubo2msa import Qubo2Msa


gurobi_available = True
try:
	import gurobipy
	import gurobi
except ImportError:
	gurobi_available = False
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


def main():
	parser = argparse.ArgumentParser("Perform MSA using D-Wave Simulator")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output_dir", required=True,
						help="The output directory")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA (will round up to nearest power of 2)")
	parser.add_argument("-l0", "--invalid_position_penalty", type=float, default=10.0,
						help="The penalty to apply for invalid sequence ordering (must be much larger than --gap_weighting)")
	parser.add_argument("-l1", "--gap_penalty", type=float, default=0.2,
						help="The penalty to apply to any detected gaps")
	parser.add_argument("-l2", "--invalid_reward_penalty", type=float, default=10.0,
						help="The penalty to apply to matches in different columns (must be greater than 1.0)")
	parser.add_argument("-d", "--delta", type=float, default=5.0,
						help="The delta factor that controls cubic to quadratic transformation (must be greater than 1.0)")
	parser.add_argument("-t", "--target", type=float, default=0.0,
						help="Allows the user to override the expected target energy of the solution (default of 0.0 means use the system defined version)")
	parser.add_argument("-n", "--repeats", type=int, default=20,
						help="This optional argument controls the argument of the same name in qbsolv.  It denotes, once a new optimal value is found, to repeat the main loop of the algorithm this number of times with no change in optimal value before stopping. The default value is 20.")
	parser.add_argument("-s", "--scoring", default=None,
						help="The scoring system to use.  Leave unset to use score of 1 for a match and 0 for mismatch.  Available options: [BLOSUM62]")
	parser.add_argument("--do_iqp", action='store_true', default=False, help="If set, run a mixed integer quadratic solver (gurobi) on the integer representation of the problem.")
	parser.add_argument("-r", "--reduced", action='store_true', default=False, help="Run in reduced mode, only E0 and E1 will be active")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	if args.do_iqp and not gurobi_available:
		print("Cannot run IQP as gurobi is not available on your system.  Please rerun without \"--do_iqp\" or install gurobi into your python environment.")

	startall = time.time()

	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)


	print("Converting problem into QUBO form")
	print("---------------------------------")
	print()

	start1 = time.time()
	m2q = Msa2Qubo(input=args.input, output=args.output_dir + "/qmsa.qubo", P=args.P, verbose=args.verbose, delta=args.delta, l0=args.invalid_position_penalty, l1=args.gap_penalty, l2=args.invalid_reward_penalty, reduced=args.reduced, scoring=args.scoring)
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
	target_energy = args.target if args.target != 0.0 else m2q.energy
	target_energy_str = "-T " + str(target_energy) if args.target else ""
	cmd_args = ['qbsolv', '-i', args.output_dir + "/qmsa.qubo", target_energy_str, "-n", str(args.repeats), "-o", args.output_dir + "/qmsa.solution"]
	#cmd_args = ['qbsolv', '-i', args.output_dir + "/qmsa.qubo", "-n", str(args.repeats), "-o", args.output_dir + "/qmsa.solution"]
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
	q2m = Qubo2Msa(settings=args.output_dir + "/qmsa.qubo.settings", solution=args.output_dir + "/qmsa.solution", input=args.input, output=args.output_dir + "/qmsa.msa", active=m2q.active, target_energy=m2q.energy, verbose=args.verbose)
	q2m.bvm(m2q.bvc)
	q2m.run()
	end3 = time.time()
	print("Time taken to interpret solution (s): ", "{0:.2f}".format(round(end3 - start3,2)))


	if args.do_iqp and gurobi_available:
		print()
		print()
		print("Solving Mixed Integer Quadratic Programming problem")
		print("---------------------------------------------------")
		print()
		gurobi.optimise(m2q.bvc)



	endall = time.time()
	print()
	print("=======================")
	print("Total time taken (s): ", "{0:.2f}".format(round(endall - startall,2)))

main()
