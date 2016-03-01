#!/usr/bin/env python3

import os
import argparse
import math
import numpy
from Bio import SeqIO


class BVC :

    __N = 0
    __Lmax = 0
    __K = 0

    __x = []

    __bvm = numpy.zeros((0, 0))

    def __init__(self, P=1, l0=0.8, l1=1.0, l2=10.0):
        self.data = []
        self.__P = P
        self.__l0 = l0
        self.__l1 = l1
        self.__l2 = l2

    def __str__(self):
        return "N=" + str(self.__N) + "\tNumber of sequences\n" \
                "Lmax=" + str(self.__Lmax) + "\tLength of longest sequence\n" \
                "K=" + str(self.__K) + "\tTotal number of characters\n" \
                "M=" + str(self.M()) + "\tLength of each row in solution space\n" \
                "m=" + str(self.m()) + "\tNumber of binary variables required to represent each position in solution space\n" \
                "P=" + str(self.__P) + "\tMaximum gap size\n" \
                "p=" + str(self.p()) + "\tNumber of binary variables required to represent each gap\n\n" \
                "l0=" + str(self.__l0) + "\tPosition weighting\n" \
                "l1=" + str(self.__l1) + "\tGap weighting\n" \
                "l2=" + str(self.__l2) + "\tReward weighting\n\n" \
                "Solution space will contain " + str(self.calc_solutionSpaceSize()) + " cells\n\n" \
                "# Binary Variables required: " + str(self.calc_minBV()) + "-" + str(self.calc_maxBV()) + "\n\n"

    def add_record(self, record):
        l = len(record.seq)
        self.__Lmax = max(self.__Lmax, l)
        self.__x.append([0] * l)      # Add a list the size of the this sequence to the 'x' list
        self.__K += l
        self.__N += 1

    def N(self):
        return self.__N

    def K(self):
        return self.__K

    def P(self):
        return self.__P

    # Returns the length of each row in solution space
    def M(self):
        return 2 * self.__Lmax + 1

    # Returns the number of binary variables required to represent each row in solution space
    def m(self):
        return math.ceil(numpy.log2(self.M()))

    # Returns the number of binary variables required to represent each gap
    def p(self):
        return math.ceil(numpy.log2(self.__P + 1))

    # Size of the solution space
    def calc_solutionSpaceSize(self):
        return self.__N * self.M()

    # Calculate lower bound on number of binary variables required
    def calc_minBV(self):
        return self.__K * (self.m() + self.p())

    # Calculate upper bound on number of binary variables required
    def calc_maxBV(self):
        return self.calc_minBV() + math.ceil(self.__N*(self.__N-1)*(self.__Lmax**2)*(1+2*self.m())/2)

    def __addE0Coefficients(self):
        x_bits = self.m()
        g_bits = self.p()
        g_offset = self.__K * self.m()
        for k in range(0, len(self.__x) - 1):
            sum_j = 0
            for j in range(0, len(self.__x[k])):
                x1_pos = k*(j+1)*x_bits
                x0_pos = k*j*x_bits
                g_pos = g_offset + (k*(j+1)*g_bits)
                for i in range(0,x_bits-1) :
                    x0i = x0_pos + i
                    x1i = x1_pos + i
                    gi = g_offset + g_pos + i
                    self.__bvm[x1i,x1i] += self.__l0
                    self.__bvm[x1i,x0i] -= self.__l0
                    #self.__bvm[x1i,?] -= self.__l0   # Linear
                    self.__bvm[x1i,gi] -= self.__l0
                    self.__bvm[x0i,x1i] -= self.__l0
                    self.__bvm[x0i,x0i] += self.__l0
                    #self.__bvm[x0i,?] -= self.__l0   # Linear
                    self.__bvm[x0i,gi] += self.__l0
                    self.__bvm[gi,x1i] -= self.__l0
                    self.__bvm[gi,x0i] += self.__l0
                    #self.__bvm[gi,?] -= self.__l0    # Linear
                    self.__bvm[gi,gi] += self.__l0

                # Coefficients for second term (gap prior to first character)
                x_pos = k*1*x_bits
                g_pos = g_offset + (k*1*g_bits)
                for i in range(0,x_bits-1) :
                    xi = x_pos + i
                    gi = g_offset + g_pos + i
                    self.__bvm[xi,xi] += self.__l0
                    self.__bvm[gi,gi] += self.__l0
                    self.__bvm[xi,gi] -= self.__l0
                    self.__bvm[gi,xi] -= self.__l0


    def __addE1Coefficients(self):
        bits = self.p()
        offset = self.__K * self.m()      # The offset required to get to the gap variables here
        for k in range(0, self.__N - 1) :
            for j in range(1, len(self.__x[k])-1) :
                g_pos = k*j*bits
                for i in range(0, bits - 1) :
                    h = offset + g_pos + i
                    self.__bvm[h,h] += self.__l1

    def createBVMatrix(self):
        size = self.calc_maxBV()
        self.__bvm = numpy.zeros((size, size))
        self.__addE0Coefficients()
        self.__addE1Coefficients()
        return self.__bvm + self.__bvm.T - numpy.diag(self.__bvm.diagonal())



# Calculates E0 ... For use in simulating the quantum computer
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
def calcE1(l1, N, G):

    sum_k = 0
    for k in range(0, N - 1) :
        sum_j = 0
        for j in range(1, len(G[k])-1):
            sum_j += G[k][j] * G[k][j]
        sum_k += sum_j

    return l1 * sum_k

# Calculates E2 ... For use in simulating the quantum computer
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


def main():

    parser = argparse.ArgumentParser("Convert Fasta file containing multiple sequences to align, into QUBO format")
    parser.add_argument("input", help="The file containing sequences to align (to be converted into QUBO)")
    parser.add_argument("-o", "--output", help="The output file, containing the QUBO representation of the MSA problem")
    parser.add_argument("-P", type=int, default=1, help="The maximum gap size allowed in the MSA")
    parser.add_argument("-l0", "--position_weighting", type=float, default=0.8, help="The weighting to apply to positioning of elements (must be larger than --gap_weighting)")
    parser.add_argument("-l1", "--gap_weighting", type=float, default=0.1, help="The weighting to apply to gap penalties")
    parser.add_argument("-l2", "--reward_weighting", type=float, default=10.0, help="The weighting to apply to reward matches (must be greater than 1.0)")
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
    bvc = BVC(args.P, l0=args.position_weighting, l1=args.gap_weighting, l2=args.reward_weighting)
    for r in records:
        bvc.add_record(r)

        # Input will only be small so just output every sequence found
        print (">" + r.id)
        print (r.seq)

    print()

    # Print current state of variables
    print(bvc)

    m = bvc.createBVMatrix()
    print (m)

main()