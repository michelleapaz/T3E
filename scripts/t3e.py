###########################################################################
# Python (v3.8.5) script t3e.py (v1.1) to run T3E algorithm
# Last update: 2022_08_25
# Author: Michelle Almeida da Paz
###########################################################################

#!/usr/bin/env python

## Import libraries
import sys
import argparse
import os
import numpy as np
import time
import random

## Set parameters
parser = argparse.ArgumentParser(description='T3E: Transposable Element Enrichment Estimator. Description: A tool for characterising the epigenetic profile of transposable elements using ChIP-seq data')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--repeat', action='store', metavar = 'repeat_file', help='transposable elements annotation [rmsk_hg38.bed (Homo sapiens) or rmsk_mm10.bed (Mus musculus)]')
parser.add_argument('--sample', action='store', metavar = 'sample_file', help='ChIP-seq sample experiment [BED format]')
parser.add_argument('--readlen', action='store', metavar = 'readlen', help='ChIP-seq input control experiment read length in base pairs [Example --readlen 36]')
parser.add_argument('--control', action='store', metavar = 'control_file', help='ChIP-seq input control experiment [BED format]')
parser.add_argument('--controlcounts', action='store', metavar = 'control_counts', help='ChIP-seq input control experiment counts [.txt format]')
parser.add_argument('--probability', action='store', metavar = 'probability_folder', help='probability folder path [Example: /control/probability/]')
parser.add_argument('--iter', action='store', metavar = 'iter', help='number of interations [Example: 100]')
parser.add_argument('--species', action='store', metavar = 'species', help='hg38 (Homo sapiens) or mm10 (Mus musculus) [Example --species hg38]')
parser.add_argument('--outputfolder', action='store', metavar = 'outputfolder', help='output folder path [Example: /probabilities]')
parser.add_argument('--outputprefix', action='store', metavar = 'outputprefix', help='prefix name of your analysis [Example: sample001]')
args = parser.parse_args()

if (len(sys.argv) == 1):
	parser.print_help()
	parser.exit()

repeats = args.repeat
sample = args.sample
read_len = int(args.readlen)
control = args.control
control_counts = args.controlcounts
prob_path = args.probability
num_iter = int(args.iter)
species = args.species
outputprefix = args.outputprefix
backup = args.outputfolder + os.path.sep + outputprefix + '_backup.txt'
background = args.outputfolder + os.path.sep + outputprefix + '_background.txt'

## Define functions
def create_chr_dict(sample):
	with open(sample, "r") as s:
		read_list = {}
		for line in s:
			line = line.rstrip()
			(chrom, start, end, read) = line.split("\t")
			read_list.setdefault(read, []).append(chrom)
	for key, value in read_list.items():
		random_chr = random.choice(value)
		if (random_chr not in chrs_dict):
			chrs_dict[random_chr] = 1
		else:
			chrs_dict[random_chr] += 1
	return chrs_dict

def read_control_file(control_file):
	array1 = []
	array2 = []
	array3 = []
	id_idx_dict = {}
	end_dict = {}
	chr_pos = ""
	first_id_seen = []
	last_seen = {}

	with open(control_file, "r") as f:
		x = 0
		value = 0
		for line in f:
			line = line.rstrip()
			(chr, pos, end, read_id) = line.split("\t")
			array2.append(pos)
			end_dict[chr] = x
			if read_id in last_seen:
				array1.append(last_seen[read_id])
				read_idx = id_idx_dict[read_id]
				array3.append(read_idx)
				map_dict[read_idx] += 1
			else:
				array1.append("-1")
				first_id_seen.append(read_id)
				id_idx_dict[read_id] = value
				read_idx = value
				value += 1
				array3.append(read_idx)
				map_dict[read_idx] = 1
			last_seen[read_id] = x
			x = x + 1
	y = 0
	for i in range(len(array1)):
		if array1[i] == "-1":
			array1[i] = last_seen[first_id_seen[y]]
			y = y + 1
		else:
			pass
	return end_dict, array1, array2, array3

def list_construction(times):
	randnum_list=[]
	rand_cumprob_list=[]
	sorted_index_prob=[]
	sorted_index_pos_list=[]
	for j in range(times):
		random = np.random.uniform(cumprob_list[0],cumprob_list[-1])
		randnum_list.append(random)
	rand_cumprob_list = list(np.append(cumprob_list, randnum_list))
	sorted_index_prob = np.argsort(rand_cumprob_list, kind="mergesort")
	sorted_index_pos_list = [index_pos_list[i] for i in sorted_index_prob]
	return sorted_index_pos_list

def rand_pos_selection():
	random_pos_list=[]
	for idx, elem in enumerate(sorted_index_pos_list):
		if (idx+1 < len(sorted_index_pos_list) and idx-1 >= 0):
			if (elem == -1):
				next_idx = idx+1

				while sorted_index_pos_list[next_idx] == -1:
					next_idx = next_idx + 1

				next_elem = sorted_index_pos_list[next_idx]
				random_pos_list.append(next_elem)
	return random_pos_list

def extract_combine_list():
	index_rand = []
	for l in range(len(random_pos_list)):
		index_rand.append(-1)
	array_comb = list(np.append(array2_2, random_pos_list))
	array_index = list(np.append(index_array2_2, index_rand))
	sorted_comb = np.argsort(array_comb, kind="mergesort")
	sorted_index = [array_index[m] for m in sorted_comb]
	return sorted_index

def access_mappings(n):
	ind_rand = 0
	list_read = []
	list_read2 = []
	read_pos = {}
	read_index = {}
	for chrom in chromosomes:
		read_pos[chrom] = []
		read_index[chrom] = []
	time_point_a = time.time()
	for i in sorted_index:
		if (i != -1):
			last = i
			list_read.append(int(last))
		else:
			random = random_pos_list[ind_rand]
			adj_random = random - int(read_len/2)

			ind_rand = ind_rand + 1
			for j in list_read:
				if int(array2_2[j]) >= int(random-read_len):
					list_read2.append(int(j))
				j += 1
			list_read = list_read2
			list_read2 = []

			if (len(list_read) == 1):
				selected = list_read[0]
			
			elif (len(list_read) > 1):
				total_conf = 0
				read_probs = []
				read_conf = []
				for elem in list_read:
					idx = int(elem) + int(start_dict[chromosome])
					read_idx = array3[idx]
					read_map = map_dict[read_idx]
					inv_map = 1/read_map
					read_conf.append(inv_map) 
					total_conf = total_conf + inv_map
				for conf in read_conf:
					read_probs.append(conf/total_conf)
				selected = np.random.choice(list_read, 1, p=read_probs)
				selected = selected[0]
			else:
				sys.exit("[ERROR] An error occurred with the list of reads")

			index = int(selected) + int(start_dict[chromosome])
			index2 = int(index)
			actual_read_pos = int(array2[index])
			diff = int(actual_read_pos) - int(adj_random)
			shift = -1 * diff
			adj_read_pos = actual_read_pos + shift
			x = 0
			while True:
				index = int(array1[index])
				for c in start_dict:
					if (index >= start_dict[c]) and (index <= end_dict[c]):
						actual_read_pos = int(array2[index])
						adj_read_pos = actual_read_pos + shift
						x = x + 1
						end_pos = int(adj_read_pos) + read_len - 1
						read_pos[c].append(int(adj_read_pos))
						read_pos[c].append(int(end_pos))
						read_index[c].append(n * (-1))
						read_index[c].append(n)
				if index == index2:
					break
			times_dict[n] = x
			n = n + 1
	return read_pos, read_index

## Main code
start_time = time.time()
print("T3E is running...")

## Global variables
if (species == "hg19" or species == "hg38"):
	chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
elif (species == "mm10"):
	chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"]
else:
	sys.exit("[ERROR] Species not defined correctly!")
start_dict={}
map_dict = {}
times_dict = {}
group_iter = {}
count_in_index = {}
group_start_chr = {}
group_end_chr = {}
group_pos = []
groups = {}
group_index = []
chrs_dict={}
chrs_dict = create_chr_dict(sample)
time_point_a = time.time()
print(f'[SAMPLE] Define the number of reads to shuffle for each chromosome: {time_point_a - start_time:.2f}s')
end_dict, array1, array2, array3 = read_control_file(control)
for elem in end_dict:
	if not start_dict:
		start_dict[elem] = 0
	else:
		start_dict[elem] = end_dict[last_chr] + 1
	last_chr = elem
time_point_b = time.time()
print(f'[CONTROL] Allocate matrix in the memory: {time_point_b - time_point_a:.2f}s')
time_point_1 = time.time()
with open(repeats, "r") as f:
	group_count = 1
	i = 0
	e = 0
	for line in f:
		line = line.replace("\n", "")
		(te_chr, te_start, te_end, repeat) = line.split("\t")
		if (te_chr not in group_start_chr):
			group_start_chr[te_chr] = i
			e = i + 1
		group_end_chr[te_chr] = e
		i = i + 2
		e = e + 2
		group_pos.append(int(te_start))
		group_pos.append(int(te_end))
		if (repeat not in groups):
			groups[repeat] = group_count
			group_index.append(group_count * (-1))
			group_index.append(group_count)
			group_count = group_count + 1
		else:
			group_index.append(groups[repeat] * (-1))
			group_index.append(groups[repeat])
time_point_2 = time.time()
print(f'[REPEATS] Read Repeat annotation: {time_point_2 - time_point_1:.2f}s')
for group in groups.values():
	count_in_index[group] = 0
for group in groups.keys():
	group_iter[group] = np.zeros(num_iter)
for chromosome in chromosomes:
	position_list = []
	cumprob_list = []
	index_list=[]
	index_pos_list=[]
	array2_2 = []
	index_array2_2 = []
	random_pos_list=[]
	sorted_index_pos_list=[]
	open_group = {}
	open_read = {}
	read_group_pos = {}
	read_group_index = {}
	sorted_read_group = {}
	sorted_read_group_pos = {}
	sorted_read_group_index = {}
	sorted_index = []
	read_pos = {}
	read_index = {}
	time_point_a = time.time()
	prob_file = prob_path + os.path.sep + chromosome + "_prob.txt"
	with open(prob_file, "r") as infile:
		for line in infile:
			line = line.rstrip()
			(pos, value) = line.split("\t")
			position_list.append(int(pos))
			cumprob_list.append(float(value))

	for t in range(chrs_dict[chromosome]):
		index_list.append(-1)
	index_pos_list = list(np.append(position_list, index_list))
	for l in array2[start_dict[chromosome]:end_dict[chromosome]+1]:
		array2_2.append(int(l))
	for l in range(len(array2_2)):
		index_array2_2.append(l)
	time_point_b = time.time()
	print(f'[PROBABILITY] Read probabilities: {time_point_b - time_point_a:.2f}s')
	for iteration in range(num_iter):
		time_sorted_index_pos_list_start = time.time()
		sorted_index_pos_list = list_construction(chrs_dict[chromosome])
		random_pos_list = rand_pos_selection()
		sorted_index = extract_combine_list()
		time_extract_combine_list_end = time.time()
		print(f'[RANDOM POSITIONS] Sample random positions: {time_extract_combine_list_end - time_sorted_index_pos_list_start:.2f}s')
		time_access_mappings_start = time.time()
		read_pos, read_index = access_mappings(group_count)
		time_access_mappings_end = time.time()
		print(f'[ACCESS MAPPINGS] Access mappings: {time_access_mappings_end - time_access_mappings_start:.2f}s')
		time_1 = time.time()
		for chromosome2 in chromosomes:
			read_group_pos[chromosome2] = []
			read_group_index[chromosome2] = []
			sorted_read_group[chromosome2] = []
			sorted_read_group_pos[chromosome2] = []
			sorted_read_group_index[chromosome2] = []
			read_group_pos[chromosome2] = list(np.append(read_pos[chromosome2], group_pos[group_start_chr[chromosome2]:group_end_chr[chromosome2]+1]))
			read_group_index[chromosome2] = list(np.append(read_index[chromosome2], group_index[group_start_chr[chromosome2]:group_end_chr[chromosome2]+1]))
			pos_index_comb = np.array(list(zip(read_group_pos[chromosome2], read_group_index[chromosome2])), dtype=[('coord', 'int32'), ('openclose', 'int32')])
			sorted_read_group[chromosome2] = np.argsort(pos_index_comb, order=['coord', 'openclose'])
			sorted_read_group_pos[chromosome2] = [read_group_pos[chromosome2][i] for i in sorted_read_group[chromosome2]]
			sorted_read_group_index[chromosome2] = [read_group_index[chromosome2][i] for i in sorted_read_group[chromosome2]]
			read_group_pos[chromosome2] = []
			read_group_index[chromosome2] = []
			sorted_read_group[chromosome2] = []
			for idx, elem in enumerate(sorted_read_group_index[chromosome2]):
				if (elem < 0):
					elem = abs(elem)
					if elem >= group_count:
						open_read.setdefault(elem, []).append(sorted_read_group_pos[chromosome2][idx])
					else:
						open_group.setdefault(elem, []).append(sorted_read_group_pos[chromosome2][idx])
				elif (elem > 0):
					elem = abs(elem)
					if (elem >= group_count):
						read_start = int(open_read[elem][0])
						read_end = sorted_read_group_pos[chromosome2][idx]
						for group in open_group.keys():
							for k_g in range(len(open_group[group])):
								group_start = int(open_group[group][k_g])
								if (open_group[group][k_g] >= read_start):	
									count_in_index[group] = count_in_index[group] + (read_end - group_start + 1) / (times_dict[elem] * read_len)
								elif (open_group[group][k_g] < read_start):
									count_in_index[group] = count_in_index[group] + (read_end - read_start + 1) / (times_dict[elem] * read_len)

						if (len(open_read[elem]) > 1):
							open_read[elem].pop(0)
						elif (len(open_read[elem]) == 1):
							open_read.pop(elem)

					elif (elem < group_count):
						group_start = int(open_group[elem][0])
						group_end = (sorted_read_group_pos[chromosome2][idx])
						for read in open_read.keys():
							for k_r in range(len(open_read[read])):
								read_start = int(open_read[read][k_r])
								if(open_group[elem][0] >= read_start):
									count_in_index[elem] = count_in_index[elem] +  (group_end - group_start + 1) / (times_dict[read] * read_len)
								elif(open_group[elem][0] < read_start):
									count_in_index[elem] = count_in_index[elem] + (group_end - read_start + 1) / (times_dict[read] * read_len)
						if (len(open_group[elem]) > 1):
							open_group[elem].pop(0)
						elif (len(open_group[elem]) == 1):
							open_group.pop(elem)
		time_2 = time.time()
		print(f'[COUNT] Count for {chromosome} - iter {iteration + 1}: {time_2 - time_1:.2f}s')
		for g_name, g_index in groups.items():
			for g_code in count_in_index:
				if (g_code == g_index):
					group_iter[g_name][iteration] = group_iter[g_name][iteration] + count_in_index[g_code]
					count_in_index[g_code] = 0
	for iteration in range(num_iter):
		with open(background, "a") as backg, open(control_counts, "r") as f:
			for line in f:
				line = line.replace("\n", "")
				(control_group_name, control_count) = line.split("\t")
				if (chromosome == "chrY"):
					print(f'iter{iteration + 1}\t{control_group_name}\t{group_iter[control_group_name][iteration]}', file=backg)
end_time = time.time()
print(end_time - start_time, " second(s)")
