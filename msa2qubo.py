#!/usr/bin/env python3

import os
import argparse
import math
import numpy
from Bio import SeqIO

# Calculates E0 ... For use in simulating the quantum computer
# E0 effectively acts as a constraint, providing high-values if invalid MSA configurations are selected
# This part allows only valid ordering by returning 0 if the configuration is valid, or signficantly greater than 0 if not
# Binary variables represent each characters position and gap size
def calcE0(l0, x, G):
    sum_k = 0
    for k in range(0, len(x) - 1):
        sum_j = 0
        for j in range(0, len(x[k])):
            term_1 = x[k][j+1] - x[k][j] - 1 - G[k][j+1]
            term_2 = x[k][0] - G[k][0]
            sum_j += (term_1 * term_1) + (term_2 * term_2)

        sum_k += sum_j

    return l0 * sum_k


# Calculates E1 ... For use in simulating the quantum computer
# E1 is used to prefer solutions that contain no gaps.  This sums all gaps across the MSA.
# The solutions with lowest numbers of gaps are prefered, but can be overcome if the rewards are great enough (see E2)
def calcE1(l1, N, G):

    sum_k = 0
    for k in range(0, N - 1) :
        sum_j = 0
        for j in range(1, len(G[k])-1):
            sum_j += G[k][j] * G[k][j]
        sum_k += sum_j

    return l1 * sum_k


# Calculates E2 ... For use in simulating the quantum computer
# E2 is used to reward alignments which match in the same column.  This function requires the a weighting matrix W is
# supplied. this matrix is pre-calculated so that characters that match will have a value of 1, those that are different
# have a value of 0.
def calcE2(l2, x, W, r):

    N = len(x)
    sum_k=0
    for k in range(0, N - 2):
        sum_q=0
        for q in range(k+1, N-1):
            sum_i=0
            for i in range(0, len(x[k])-1):
                sum_j=0
                for j in range(0, len(x[q])-1):
                    sum_j += W[i][j][k][q] * r[i][j][k][q] * math.floor(1 - l2 * math.pow(x[k][i] - x[q][j], 2))
                sum_i += sum_j
            sum_q += sum_i
        sum_k += sum_q

    return -sum_k


def createW(x, Lmax):

	N=len(x)
	W=numpy.empty(shape=(Lmax,Lmax,N,N))

	for k in range(0,N-2):
		for q in range(k+1,N-1):
			for i in range(0,len(x[k])-1):
				for j in range(0,len(x[q])-1):
					W[i][j][k][q] = 1 if x[k][i] == x[q][j] else 0

	return W

def countsW(W, x):

	z=0
	nz=0
	N=len(x)
	for k in range(0,N-2):
		for q in range(k+1,N-1):
			for i in range(0,len(x[k])-1):
				for j in range(0,len(x[q])-1):
					if W[i][j][k][q] > 0:
						nz += 1
					else:
						z+=1

	return nz, z


def main():

	parser = argparse.ArgumentParser("Convert Fasta file containing multiple sequences to align, into QUBO format")
	parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
	parser.add_argument("-o", "--output", help="The output file, containing the QUBO representation of the MSA problem")
	parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA")
	parser.add_argument("-l0", "--position_weighting", type=float, default=0.8, help="The weighting to apply to positioning of elements (must be larger than --gap_weighting)")
	parser.add_argument("-l1", "--gap_weighting", type=float, default=0.1, help="The weighting to apply to gap penalties")
	parser.add_argument("-l2", "--reward_weighting", type=float, default=10.0, help="The weighting to apply to reward matches (must be greater than 1.0)")
	parser.add_argument("-d", "--delta", type=float, default="2.0", help="The weight to apply when converting cubic reward terms to quadratic")
	parser.add_argument("-v", "--verbose", help="Display extra information")
	args = parser.parse_args()

	if args.position_weighting < args.gap_weighting:
	    print("-l0 must be larger than -l1")
	    exit(1)

	if args.gap_weighting <= 0.0:
	    print("-l1 must be larger than 0.0")
	    print(1)

	if args.reward_weighting <= 1.0:
	    print("-l2 must be > 1.0")
	    exit(1)


	print("Loading input into memory...", end="")
	handle = open(args.input, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	print(" done")
	print()
	print("Input:")

	# Collect variables
	N = len(records)    # Number of sequences
	Lmax = 0            # Length of longest sequence
	x = []              # Solution variables defining Base10 position of each character in the MSA (2D list)
	c = []
	K = 0               # Total number of characters
	for r in records:
	    l = len(r.seq)
	    c.append(r.seq)
	    Lmax = max(Lmax, l)
	    x.append([0] * l)      # Add a list the size of the this sequence to the 'x' list
	    K += l

	    # Input will only be small so just output every sequence found
	    print (">" + r.id)
	    print (r.seq)

	M = 2*Lmax+1        # Maximum length of a row in solution space
	m = math.ceil(numpy.log2(M))   # Number of binary variables required to represent each position

	p = math.ceil(numpy.log2(args.P + 1))  # Number of binary variables required to represent each gap

	min_bv = K * (m + p)
	max_bv = min_bv + math.ceil(N*(N-1)*(Lmax**2)*(1+2*m)/2)

	print()
	print("Variables:")
	print("N=" + str(N) + "\t- Number of sequences")
	print("Lmax=" + str(Lmax) + "\t- Length of longest sequence")
	print("K=" + str(K) + "\t- Total number of characters in MSA")
	print("M=" + str(M) + "\t- Maximum length of a row in solution space")
	print("m=" + str(m) + "\t- Number of binary variables required to represent the position of each character")
	print("P=" + str(args.P) + "\t- Maximum gap size")
	print("p=" + str(p) + "\t- Number of binary variables required to represet each gap")
	print()
	print("l0=" + str(args.position_weighting) + "\t- Position weighting")
	print("l1=" + str(args.gap_weighting) + "\t- Gap weighting")
	print("l2=" + str(args.reward_weighting) + "\t- Reward weighting")
	print("d=" + str(args.delta) + "\t- Cubic to quadratic weighting")
	print()
	print("Solution space will be performed on a grid of " + str(N) + " rows and " + str(M) + " columns")
	print()
	print("Estimates for total number of binary variables required (TBV). " + str(min_bv) + " <= TBV <= " + str(max_bv))
	print()

	# Binary variables for E0 (positioning options)
	#b = [K*m]
	#g = [K*p]
	#E0_bv = []
	#E0_bv.extend(b)
	#E0.bv.extend(g)
	print("Binary variables required for x: " + str(K*m))
	print("Binary variables required for G: " + str(K*p))
	print()

	W=createW(c, Lmax)

	print("W contains " + str(W.size) + " elements of shape " + str(W.shape))
	print()

	nz,z=countsW(W, x)
	print("However, not all elements are used in W.  W contains " + str(nz) + " non-zero elements (to be used as cubic binary variables) and " + str(z) + " elements, totalling " + str(nz+z) + " considered elements")

	qbv = nz*2
	print("After substituting to quadratic form, the E2 part of the energy function requires " + str(qbv) + " binary variables")
	print()
	print("Total number of binary variables required for input: " + str(qbv+min_bv))


main()