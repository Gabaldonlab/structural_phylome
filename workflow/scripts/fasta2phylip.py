#!/usr/bin/env python
# adapted from https://github.com/jvollme/fasta2phylip/blob/master/fasta2phylip.py
import os
from Bio import AlignIO

def main():
	infile = open(snakemake.input[0], "r")
	outfile = open(snakemake.output[0], "w")
	alignments = AlignIO.parse(infile, "fasta")
	AlignIO.write(alignments, outfile, "phylip-relaxed")
	infile.close()
	outfile.close()
	
main()
