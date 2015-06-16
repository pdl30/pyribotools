#!/usr/bin/python

########################################################################
# 27 April 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools
import HTSeq
from multiprocessing import Pool, Manager
import argparse

def find_index(index_no):
	with open("/home/patrick/Scripts/pyribotools/adapter_indices_table.txt") as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if index_no == int(word[1]):
				return word[3]

def cut(fq, conditions, outdir, index, rpi):
	name = os.path.basename(fq)
	index = find_index(int(conditions[fq]))
	command = "cutadapt -O 3 -b {} -a {} {} -q 20 -m 5 > {}/{}".format(rpi, index, fq, outdir, name)
	subprocess.call(command, shell=True)

def strip_first_base(fq, idir, outdir):
	name = os.path.basename(fq)
	ifile = idir + '/' + name
	name = re.sub(".fq", ".fa", name)
	output = open(outdir + '/' + name, "w")
	for s in HTSeq.FastqReader( ifile ):
		#output.write(">{}\n{}\n+\n{}\n".format(s.name, s.seq[1:], s.qualstr[1:])),
		output.write(">{}\n{}\n".format(s.name, s.seq[1:])),
	output.close()

def remove_rrna(fq, idir, outdir, rrna_reference):
	name = os.path.basename(fq)
	name = re.sub(".fq", ".fa", name)
	ifile = idir + '/' + name
	output = outdir + '/' + name
	command = "bowtie -f --seedlen=23 --un={} {} {} >/dev/null".format(output, rrna_reference, ifile)
	subprocess.call(command, shell=True)

def run_alignment(fq, idir, outdir, gtf, index):
	name = os.path.basename(fq)
	name = re.sub(".fq", ".fa", name)
	ifile = idir + '/' + name
	sample_name=  re.sub(".fa", "", name)
	output = outdir + '/' + sample_name
	command = "tophat --no-novel-juncs -p 6 --output-dir {} --GTF {} {} {}".format(output, gtf, index, ifile)
	subprocess.call(command.split())

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

def f1(args):
	return cut(*args)
def f2(args):
	return strip_first_base(*args)
def f3(args):
	return remove_rrna(*args)
def f4(args):
	return run_alignment(*args)

def main():
	parser = argparse.ArgumentParser(description='Will create trimmed_fqs, first_char_stripped and alignment_dir\n')
	parser.add_argument('-c', '--config', help='Contains [Conditions] with Fastq and index as value', required=True) 
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/hg19', required=True)
	parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False) 
	parser.add_argument('-o', '--output', help='Output place for directories', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	

	rpi = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA";

	if args["genome"] == "hg19":
		rrna_reference = "/home/patrick/Reference_Genomes/rRNA_indices/human/rrna_index/rrna_reference"
		gtf = "/home/patrick/Reference_Genomes/pyngspipe_references/hg19/hg19.gtf"
		index = "/home/patrick/Reference_Genomes/pyngspipe_references/hg19/hg19"
	else:
		rrna_reference = "/home/patrick/Reference_Genomes/rRNA_indices/mouse/rrna_index/rrna_index_mouse"
		gtf = "/home/patrick/Reference_Genomes/mm10/Ensembl/76/Mus_musculus.GRCm38.76_ucsc.gtf"
		index = "/home/patrick/Reference_Genomes/pyngspipe_references/mm10/mm10"
	
	trimmed_dir = args["output"] + '/trimmed_fqs'
	first_char_stripped_dir = args["output"]+ '/first_char_stripped'
	norna_dir = args["output"]+ '/norrna'
	align_dir = args["output"]+ '/alignment'
	if os.path.isdir(trimmed_dir):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", trimmed_dir])
	if os.path.isdir(first_char_stripped_dir):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", first_char_stripped_dir])
	if os.path.isdir(norna_dir):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", norna_dir])
	if os.path.isdir(align_dir):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", align_dir])

	pool = Pool(int(args["threads"]))
	pool.map(f1, itertools.izip(list(conditions.keys()), itertools.repeat(conditions), itertools.repeat(trimmed_dir), itertools.repeat(index), itertools.repeat(rpi)))
	pool.map(f2, itertools.izip(list(conditions.keys()), itertools.repeat(trimmed_dir), itertools.repeat(first_char_stripped_dir)))
	pool.map(f3, itertools.izip(list(conditions.keys()), itertools.repeat(first_char_stripped_dir), itertools.repeat(norna_dir), itertools.repeat(rrna_reference)))
	pool.map(f4, itertools.izip(list(conditions.keys()), itertools.repeat(norna_dir), itertools.repeat(align_dir), itertools.repeat(gtf), itertools.repeat(index)))
