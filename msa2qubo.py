#!/usr/bin/env python3

"""msa2qubo.py: Converts a fasta file representing a multiple sequence alignment into a QUBO format file."""

import argparse
import sys
from Bio import SeqIO
from bvc import BVC
import weightings

__author__ = "Dan Mapleson, Luis Yanes, Katie Barr, Sophie Kirkwood and Tim Stitt"
__copyright__ = "Copyright 2016, Quantum MSA"
__credits__ = ["Dan Mapleson", "Luis Yanes", "Katie Barr",
                    "Sophie Kirkwood", "Tim Stitt"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Dan Mapleson,"
__email__ = "daniel.mapleson@earlham.ac.uk"
__status__ = "Prototype"


class Msa2Qubo:

	def __init__(self, input, output, P, delta, l0, l1, l2, reduced=False, verbose=False, scoring = None):
		self.data = []
		self.input = input
		self.output = output
		self.P = P
		self.delta = delta
		self.l0 = l0
		self.l1 = l1
		self.l2 = l2
		self.reduced = reduced
		self.verbose = verbose
		self.scoring = scoring

		self.energy = 0.0
		self.active = []
		self.bvc = None


	def run(self):

		if self.l0 < self.l1:
			print("-l0 must be larger than -l1")
			exit(1)

		if self.l1 <= 0.0:
			print("-l1 must be larger than 0.0")
			print(1)

		if self.l2 <= 1.0:
			print("-l2 must be > 1.0")
			exit(1)

		if self.delta <= 1.0:
			print("-d must be > 1.0")
			exit(1)

		if self.P < 1:
			print("-P must be >= 1")
			exit(1)

		weights = None
		if self.scoring:
			print("Loading " + self.scoring + " scoring matrix...", end="")
			weights = weightings.loadScoring(self.scoring)
			print(" done")


		print("Loading input into memory...", end="")
		handle = open(self.input, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		print(" done")
		print()
		print("Input:")

		# Collect variables
		self.bvc = BVC(P=self.P, d=self.delta, l0=self.l0, l1=self.l1, l2=self.l2, reduced=self.reduced, weights=weights)
		for r in records:
			self.bvc.add_record(r)

			# Input will only be small so just output every sequence found
			print(">" + r.id)
			print(r.seq)

		print()

		# Print current state of variables
		if self.verbose:
			print(self.bvc)

		if self.bvc.p() >= self.bvc.m():
			print("P (max gap size) is too large.  Please ensure that p (bits per gap) is significantly less than m (bits per position in solution space)")
			exit(1)

		# Create matrix
		print()
		print("Creating W matrix ...", end="")
		sys.stdout.flush()
		self.bvc.createW()
		# print (m)
		print(" done")

		self.bvc.createBVMatrix(intmode=True, verbose=self.verbose)

		if self.verbose:
			self.bvc.printIntegerCoefficients()

		print("Creating matrix of coefficients of binary variables ...", end="")
		sys.stdout.flush()
		self.bvc.createBVMatrix(verbose=self.verbose)
		print(" done")
		print("Number of active binary variables:", self.bvc.calcActiveBVs())
		self.active = self.bvc.active

		self.energy = -self.bvc.energy
		print("Target energy:", self.energy)

		# Write QUBO file to disk
		print("Writing QUBO output to disk ...", end="")
		sys.stdout.flush()
		self.bvc.writeQUBO(self.output, self.input)
		print(" done")
		print("Output saved to:" + self.output)
		sys.stdout.flush()

		if self.verbose:
			self.bvc.plotMatrix(self.output + ".png")

		# Save settings file
		print()
		self.bvc.save_settings(self.output + ".settings")
		print("Saved settings to: " + self.output + ".settings")



def main():
	parser = argparse.ArgumentParser("Convert Fasta file containing multiple sequences to align, into QUBO format")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output", required=True,
						help="The output file, containing the QUBO representation of the MSA problem")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA (will round up to nearest power of 2)")
	parser.add_argument("-d", "--delta", type=float, default=2.0, help="Delta.  The scaling factor to apply when converting products of 3BVs to 2BVs.")
	parser.add_argument("-l0", "--invalid_position_penalty", type=float, default=10.0,
						help="The penalty to apply for invalid sequence ordering (must be much larger than --gap_weighting)")
	parser.add_argument("-l1", "--gap_penalty", type=float, default=0.1,
						help="The penalty to apply to any detected gaps")
	parser.add_argument("-l2", "--invalid_reward_penalty", type=float, default=100.0,
						help="The penalty to apply to matches in different columns (must be greater than 1.0)")
	parser.add_argument("-r", "--reduced", action='store_true', default=False,
						help="Run in reduced mode, only E0 and E1 active")
	parser.add_argument("-s", "--scoring", default=None,
						help="The scoring system to use.  Leave unset to use score of 1 for a match and 0 for mismatch.  Available options: [BLOSUM62]")
	parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Display extra information")
	args = parser.parse_args()

	if args.invalid_position_penalty < args.gap_penalty:
		print("-l0 must be larger than -l1")
		exit(1)

	if args.gap_penalty <= 0.0:
		print("-l1 must be larger than 0.0")
		print(1)

	if args.invalid_reward_penalty <= 1.0:
		print("-l2 must be > 1.0")
		exit(1)

	if args.delta <= 1.0:
		print("-d must be > 1.0")
		exit(1)

	if args.P < 1:
		print("-P must be >= 1")
		exit(1)

	m2q = Msa2Qubo(input=args.input, output=args.output, P=args.P, delta=args.delta, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting,
			reduced=args.reduced, scoring=args.scoring, verbose=args.verbose)
	m2q.run()


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
