###########################################################################
# Python (v3.8.5) script enrichment.py (v1.0) to calculate enrichments
# Last update: 2022_06_07
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/env python

import math
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='T3E: Enrichment')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--background', action='store', metavar = 'background', help='background file created by T3E [Example: sample001_background.txt]')
parser.add_argument('--signal', action='store', metavar = 'signal', help='ChIP-seq sample experiment counts [.txt format]')
parser.add_argument('--iter', action='store', metavar = 'iter', help='number of interations [Example: 100]')
parser.add_argument('--alpha', action='store', metavar = 'alpha', help='level of significance to report enrichment [Example: 0.05]')
parser.add_argument('--enrichment', action='store', metavar = 'enrichment', help='log2FC threshold to report enrichment [Example: 1.0]')
parser.add_argument('--outputfolder', action='store', metavar = 'outputfolder', help='output folder path [Example: /results]')
parser.add_argument('--outputprefix', action='store', metavar = 'outputprefix', help='prefix name of your analysis [Example: sample001]')
args = parser.parse_args()

if (len(sys.argv) == 1):
	parser.print_help()
	parser.exit()

background = args.background
signal = args.signal
num_iter = int(args.iter)
alpha = float(args.alpha)
enrichment = float(args.enrichment)
outputprefix = args.outputprefix
output = args.outputfolder + os.path.sep + outputprefix + '_enrichment.txt'

def open_signal(sample):
	with open(sample, "r") as s:
		for line in s:
			line = line.rstrip()
			(repeat, counts) = line.split("\t")
			repeats_counts[repeat] = counts
			qt[repeat] = 0
			sum_backg[repeat] = 0
	return repeats_counts
	
def open_iterations(iterations):
	with open(iterations, "r") as i:
		for line in i:
			line = line.rstrip()
			(iterate, repeat, counts) = line.split("\t")
			if repeat in repeats_counts:
				if (float(counts) > float(repeats_counts[repeat])):
					qt[repeat] += 1
				sum_backg[repeat] += float(counts)
	return qt, sum_backg

repeats_counts = {}
qt = {}
sum_backg = {}
repeats_counts = open_signal(signal)

qt, sum_backg = open_iterations(background)

with open(output, "w") as o:
	for repeat in sum_backg.keys():
		pvalue = qt[repeat]/num_iter
		mean = sum_backg[repeat]/num_iter
		if (mean > 0):
			foldchange = float(repeats_counts[repeat])/mean
			log2fc = math.log(foldchange,2)
			if ((pvalue <= alpha) and (log2fc >= enrichment)):
				print(repeat, pvalue, log2fc)
			print(repeat, "\t", pvalue, "\t", log2fc, file=o)
