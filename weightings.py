#!/usr/bin/env python3

"""weightings.py: Contains some commonly used scoring matrices that can be used to weight alignments"""

import os
import sys

def loadScoring(scoringtype):

	here = os.path.dirname(os.path.realpath(__file__))

	scoring_file = os.path.join(here, "scoring", scoringtype + ".txt")

	matrix = {}
	keys = []

	with open(scoring_file) as f:
		header = True
		for line in f:

			line = line.strip()

			if line and not line.startswith("#"):	# If c is not an empty line and isn't a comment then we have something to work with

				parts = line.split()

				if header:
					for p in parts:
						matrix[p] = None
						keys.append(p)
					header = False
				else:
					m = {}
					for i in range(len(parts)-1):
						k = keys[i]
						m[k] = float(parts[i+1])
					matrix[parts[0]] = m

		maxval = -10000000
		minval = 10000000
		for ko, vo in matrix.items():
			for ki, vi in matrix[ko].items():
				maxval = max(maxval, vi)
				minval = min(minval, vi)

		diff = maxval - minval

		for ko, vo in matrix.items():
			for ki, vi in matrix[ko].items():
				matrix[ko][ki] = float(vi - minval) / float(diff)
	return matrix



