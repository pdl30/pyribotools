#!/usr/bin/python

########################################################################
# 15 June 2015
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

def align(conditions, index, outdir):
	for fq in conditions:
		name = os.path.basename(fq)
		name = re.sub(".fa", "", name)
		f = open("{}/{}_report.txt".format(outdir, name), "w")
		command = "bowtie -f -m 1 -v 1 --al {3}/{1}.fa -p 20 -a --best --strata --sam {0} {2} > {3}/{1}.sam".format(index, name, fq, outdir)
		subprocess.call(command, shell=True, stderr=f)
		command2 = "bowtie -m 500 -v 3 -p 20 -a --best --strat -f --sam {0} {1}/{2}.fa > {1}/{2}_utr.sam".format(index, outdir, name)
		subprocess.call(command2, shell=True, stderr=f)

def convert_bam(conditions, outdir):
	for fq in conditions:
		name = os.path.basename(fq)
		name = re.sub(".fa", "", name)
		command = "samtools view -bS {0}/{1}_utr.sam > {0}/{1}_utr.bam".format(outdir, name)
		subprocess.call(command, shell=True)
		command = "samtools sort {0}/{1}_utr.bam {0}/{1}_utr_sort".format(outdir, name)
		subprocess.call(command, shell=True)
		command = "samtools index {}/{}_utr_sort.bam".format(outdir, name)
		subprocess.call(command, shell=True)

def read_bam(fq, outdir, return_dict, genes): 
	profile = numpy.zeros( 2*200, dtype="f" )
	constant = 1
	name = os.path.basename(fq)
	name = re.sub(".fa", "", name)
	bamfile = HTSeq.BAM_Reader("{}/{}_utr_sort.bam".format(outdir, name))
	for almnt in bamfile:
		if almnt.iv:
			if almnt.iv.strand == "+":
				start_in_window = almnt.iv.start# - mid + halfwinwidth
				end_in_window   = almnt.iv.end   #- mid +  halfwinwidth
				start_in_window = max( start_in_window, 0 ) #A more elegant solution for bed files under 0
				end_in_window = min( end_in_window, 400 )
				if start_in_window >= 400 or end_in_window < 0:
					continue
				profile[ start_in_window : end_in_window ] += constant
				if (almnt.iv.chrom, fq) not in genes:
					genes[(almnt.iv.chrom, fq)] = 1
	return_dict[fq] = profile

def function1(args):
	return read_bam(*args)

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

def plot_utrs(conditions, rev_conds, outdir):
	pool = Pool(24)
	manager = Manager()
	return_dict = manager.dict()
	genes = manager.dict()
	pool.map(function1, itertools.izip(list(conditions.keys()), itertools.repeat(outdir), itertools.repeat(return_dict), itertools.repeat(genes)))
	combined_profiles = {}
	normal = {}
	#Have to add all conditions together
	today = date.today()
	date_format = "{}_{}_{}".format(today.day, today.month, today.year)
	pp = PdfPages("{}/{}_UTR_averaged.pdf".format(outdir, date_format))

	#First need to average over all genes:
	averaged_profiles = {}
	len_genes = {}
	for key1, key2 in genes.keys(): #Chromosome, bam
		if key2 not in len_genes:
			len_genes[key2] = {}
			len_genes[key2][key1] = 1
		else:
			len_genes[key2][key1] = 1

	for key in return_dict.keys():
		averaged_profiles[key] = return_dict[key]/len(len_genes[key].keys())
	
	for key in rev_conds:
		fig = pyplot.figure()
		pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green"])
		for fasta in rev_conds[key]:
			name = re.sub(".fa", "", fasta)
			name = os.path.basename(name)
			uniq_count = read_reports('{}/{}_report.txt'.format(outdir, name))
			normal[fasta] = uniq_count
			norm = 100000/float(uniq_count)
			normalised_profile = norm * averaged_profiles[fasta]
			if key not in combined_profiles:
				combined_profiles[fasta] = normalised_profile
			else:
				combined_profiles[fasta] += normalised_profile

			pyplot.plot( numpy.arange( -200, 200 ), normalised_profile, label=name)
		pyplot.legend(prop={'size':6})
		pyplot.title(key)  
		pp.savefig(fig)
		pyplot.close()
	fig = pyplot.figure()
	pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green"])
	for key in combined_profiles:
		name = re.sub(".fa", "", key)
		name = os.path.basename(name)
		pyplot.plot( numpy.arange( -200, 200 ), combined_profiles[key], label=name)  
	pyplot.legend(prop={'size':6})
	pp.savefig(fig)
	pp.close()

def ConfigSectionMap(Config, section):
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

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def get_config_args(args, arglist):
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	dict_list = []
	for arg in arglist:
		dict_list.append(ConfigSectionMap(Config, arg))
	return dict_list

def main():
	parser = argparse.ArgumentParser(description='UTR profiling for Riboseq Analysis.\n')
	parser.add_argument('-c','--config', help='Config file containing [Conditions], these are fastqs and conditions', required=True)
	parser.add_argument('-i','--index', help='Bowtie UTR index', required=True)
	parser.add_argument('-o', '--outdir', help='Output directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	conditions = get_config_args(args, ["Conditions"])[0]
	rev_conds = reverse_dict(conditions)
	if not os.path.isdir(args["outdir"]):
		os.mkdir(args["outdir"])
	align(conditions, args["index"], args["outdir"])
	convert_bam(conditions, args["outdir"])
	rev_conds = reverse_dict(conditions)
	plot_utrs(conditions, rev_conds, args["outdir"])
