#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import argparse
import HTSeq
import numpy
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool, Manager
import itertools
import ConfigParser
import tempfile
from datetime import date
from mpltools import color

def read_bam_counts(bam, counts_dict):
	bamfile = HTSeq.BAM_Reader(bam)
	constant = 1
	name = re.sub("_sort.bam", "", bam)
	for almnt in bamfile:
		if almnt.iv:
			if almnt.iv.strand == "+":
				start_in_window = almnt.iv.start# - mid + halfwinwidth
				end_in_window   = almnt.iv.end   #- mid +  halfwinwidth
				start_in_window = max( start_in_window, 0 ) #A more elegant solution for bed files under 0
				end_in_window = min( end_in_window, 400 )
				if start_in_window >= 400 or end_in_window < 0:
					continue
				if (almnt.iv.chrom, bam) not in counts_dict:
					counts_dict[(almnt.iv.chrom, bam)] = 1
				else:
					counts_dict[(almnt.iv.chrom, bam)] += 1

def read_reports(report):
	with open(report) as f:
		i = 0
		for line in f:
			if i == 6:
				line=line.rstrip()
				line = line.strip("# reads with at least one reported alignment: ")
				word = line.split()
				return int(word[0])
			i += 1

def get_counts(conditions, rev_conds, outdir):
	pool = Pool(24)
	manager = Manager()
	counts_dict = manager.dict()
	split_counts_dict = manager.dict()
	pool.map(function2, itertools.izip(list(conditions.keys()), itertools.repeat(counts_dict), itertools.repeat(split_counts_dict)))
	output = open("{}/utr_counts.tsv".format(human_count_results), "w")
	for cond in sorted(conditions):
		name = os.path.basename(cond)
		output.write("\t{}".format(name)),
	output.write("\n"),
	normal = {}
	#FORMAT THE OUTPUT DICTS
	counts_dict2 = {}
	for key1, key2 in counts_dict.keys(): #Key1 == chromosome key2 == bam file
		if key1 not in counts_dict2:
			counts_dict2[key1] = {}
			counts_dict2[key1][key2] = counts_dict[(key1, key2)]
		else:
			counts_dict2[key1][key2] = counts_dict[(key1, key2)]
	#Write normal UTR output
	for utr in counts_dict2.keys(): 
		output.write("{}".format(utr)),
		for bam in sorted(conditions):
			name = re.sub("_utr_sort.bam", "", bam)
			name2 = os.path.basename(name)
			uniq_count = read_reports(name + '_report.txt')
			normal[bam] = uniq_count
			norm = 100000/float(normal[bam])
			count = counts_dict2[utr].get(bam, 0)
			if count > 0:
				norm2 = norm * count
			else:
				norm2 = 0
			output.write("\t{}".format(norm2)),
		output.write("\n"),
	output.close()

def function2(args):
	return read_bam_counts(*args)

def create_design_for_R(idict):
	fh = tempfile.NamedTemporaryFile(delete = False)
	fh.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		bam_name = os.path.basename(key)
		name = re.sub("_sort.bam$", "", bam_name)
		fh.write("{}\t{}\t{}\n".format(name, bam_name, idict[key]))
	fh.close()
	return fh.name

def run_deseq(conditions, rev_conds, outdir):
	today = date.today()
	date_format = "{}_{}_{}".format(today.day, today.month, today.year)
	pdata = create_design_for_R(conditions)
	combinations = list(itertools.combinations(rev_conds.keys(),2))
	for key in combinations:
		command = "Rscript {}/deseq2_rcode2.R {}/human_utr_counts.tsv {} {} {} {}/human_plots.pdf {}".format(scripts_dir, human_count_results, pdata, key[0], key[1], human_utr_deseq, human_utr_deseq)
		subprocess.Popen(command.split())

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def get_config_args(args, arglist):
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	dict_list = []
	for arg in arglist:
		dict_list.append(ConfigSectionMap(Config, arg))
	return dict_list

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def main():
	parser = argparse.ArgumentParser(description='DESEQ UTR Analysis.\n')
	parser.add_argument('-c','--config', help='Config file containing [Conditions]', required=True)
	parser.add_argument('-o', '--outdir', help='Output directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	conditions = get_config_args(args, ["Conditions"])[0]
	rev_conds = reverse_dict(conditions)
	if not os.path.isdir(args["outdir"]):
		os.mkdir(args["outdir"])
	get_counts(conditions, rev_conds, args["outdir"])
	run_deseq(conditions, rev_conds, args["outdir"])