#MSA via QUBO

These scripts convert a Fasta format input file representing sequences to align against each other.  The qmsa.py script
first calls msa2qubo.py which converts this input into a QUBO representation of a multiple sequence alignment algorithm.  The QUBO file generated can
then be run on a D-Wave simulator or on a D-Wave device itself.  The results of the D-Wave (simulator) is a list of
binary numbers, which are then interpreted using qubo2msa.py to show the results of the sequence alignment.

##Installation:

  - Install python3 with numpy (python_anaconda is recommended)
  - Add biopython to your python installation
  - (Optional) Install qOp and add the bin directory to your PATH.
  - (Optional) Install gurobi for integer quadratic programming option.


##Operating Instructions:

Type ```./qmsa.py --help``` for instructions on how to run the script.  The user is presented with options to specify lambda0,
lamdba1 and lambda2.

##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

 * Daniel Mapleson
 * Luis Yanes
 * Katie Barr
 * Sophie Kirkwood
 * Tim Stitt


##Acknowledgements:

D-Wave
NASA
