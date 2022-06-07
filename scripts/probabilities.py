###########################################################################
# Python (v3.8.5) script probabilities.py (v1.0) for input-based probabilities
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/env python

import sys
import argparse
import os
import time

parser = argparse.ArgumentParser(description='Calculate input-basd background probability distribution')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--control', action='store', metavar = 'control_file', help='ChIP-seq input control experiment [BED format]')
parser.add_argument('--readlen', action='store', metavar = 'readlen', help='ChIP-seq input control experiment read length in base pairs [Example --readlen 36]')
parser.add_argument('--species', action='store', metavar = 'species', help='hg38 (Homo sapiens) or mm10 (Mus musculus) [Example --species hg38]')
parser.add_argument('--outputfolder', action='store', metavar = 'outputfolder', help='output folder path [Example: /probabilities]')
args = parser.parse_args()

if (len(sys.argv) == 1):
	parser.print_help()
	parser.exit()

control = args.control
read_len = int(args.readlen)
species = args.species
outputfolder = args.outputfolder

def read_frequency():
	time_point_a = time.time()
	previous_start = 0
	previous_end = 0
	with open(control, "r") as f:
		for line in f:
			line = line.replace("\n", "")
			(chrom, start, end, read) = line.split("\t")
			if (read not in freq_reads):
				freq_reads[read] = 1
			else:
				freq_reads[read] += 1
			chr_reads[chrom].append(read)
	time_point_b = time.time()
	print(f'#Time for read freq: {time_point_b - time_point_a:.5f}s')

	return freq_reads

def calculate_probability(c):
	time_point_a = time.time()
	with open(control, "r") as f:
		for line in f:
			line = line.replace("\n", "")
			(chrom, start, end, read) = line.split("\t")
			if (chrom == c):
				for i in range(int(start), int(end)):
					if (i not in pos_prob):
						pos_prob[i] = (1/freq_reads[read])/(read_len*num_reads[chrom])
					else:
						pos_prob[i] += (1/freq_reads[read])/(read_len*num_reads[chrom])
	time_point_b = time.time()
	print(f'#Time for calculate probability: {time_point_b - time_point_a:.5f}s')

	return pos_prob

if (species == "hg19" or species == "hg38"):
	chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
elif (species == "mm10"):
	chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"]
else:
	sys.exit("[ERROR] Species not defined correctly!")

pos_prob = {}
num_reads = {}
chr_reads = {}
freq_reads = {}

for chromosome in chromosomes:
	num_reads[chromosome] = 0
	chr_reads[chromosome] = []

freq_reads= read_frequency()

for chromosome in chromosomes:
	for elem in chr_reads[chromosome]:
		num_reads[chromosome] = num_reads[chromosome] + 1/freq_reads[elem]
	pos_prob = calculate_probability(chromosome)
	first_pos = int(list(pos_prob.keys())[0])
	file_out = args.outputfolder + os.path.sep + chromosome + '_prob.txt'
	with open(file_out, "w") as o:
		for pos in pos_prob:
			if(pos == first_pos):
				print(pos, "\t", pos_prob[pos], file=o)
				last_cumprob = pos_prob[pos]
			else:
				cumprob = pos_prob[pos] + last_cumprob
				print(pos, "\t", cumprob, file=o)
				last_cumprob = cumprob
	pos_prob = {}
