# -*- coding: utf-8 -*-
import regex as re
import pandas as pd
import argparse, os, logging, sys, traceback, time, glob, socket
import numpy as np
from Bio import Entrez, SeqIO, Seq
import PA_limitedIndel as PA_script
from Bio.Alphabet import generic_dna
from subprocess import Popen, PIPE
import CRISTA_main

GENOMES_PATH = "/groups/itay_mayrose/shiranabad/CRISPR/bwa_genomes/"

N_GAPS = 0#3
N_DISTEDIT = 4
CHROMOSOME_HEADING = "chromosome"
STARTPOS_HEADING = "start position"
ENDPOS_HEADING = "end position"
STRAND_HEADING = "strand"
SGRNA_HEADING = "aligned sgRNA"
OFFTARGET_HEADING = "aligned site"
MATCH_SCORE = 1
MISMATCH_PENALTY = 0
GAP_PENALTY = -1.25
EXTENSION = 100
#SAMTOOLS_EXE = "/groups/itay_mayrose/shiranabad/applications/samtools-1.2/samtools "
SAMTOOLS_EXE = "/groups/itay_mayrose/shiranabad/applications/samtools-1.9/samtools " # updated for power server
#HOSTNAME = socket.gethostname()
# if HOSTNAME == 'jekyl.tau.ac.il' or HOSTNAME.startswith("compute-7-") or HOSTNAME.startswith("compute-8-"):
# 	BWA_EXE = BWA_COMMAND_ON_SERVER = "/share/apps/bin/bwa"
# else:
# 	BWA_EXE = BWA_COMMAND_ON_SERVER = "/share/apps/bwa-0.7.15/bwa"
BWA_EXE = BWA_COMMAND_ON_SERVER = "bwa"


def get_seq_by_orientation(seq, strand):
	"""
		returns the original sequence if it's the plus strand, o/w returns the reverse-complement
	"""
	if strand == "-":
		seqobject = Seq.Seq(seq, generic_dna)
		seq = str(Seq.Seq.reverse_complement(seqobject))
	return seq


def fetch_dna_coordinates_offline(genome_assembly, chrom, startpos, endpos, target_site=None):
	"""
	reads the sequence from the local copy of the genome assembly. assumes that the sequence does not spread over more than 2 chop files
	:return:
	"""
	seq = fetch_dna_coordinates_faidx(genome_assembly, chrom=chrom[3:], start=startpos, end=endpos)
	if seq == "":
		seq = fetch_dna_coordinates_faidx(genome_assembly, chrom=chrom, start=startpos, end=endpos)
	return seq


def fetch_dna_coordinates_faidx(genome_assembly, chrom, start, end):
	"""Call tabix and generate an array of strings for each line it returns."""
	filename = get_bwa_genome_path(genome_assembly)
	if not os.path.exists(filename):
		raise ValueError("Genome assembly is not valid. Please contact us for assistance.")
	query = '{}:{}-{}'.format(chrom, start, end)
	#print(SAMTOOLS_EXE + ' faidx ', filename, " ", query)
	process = Popen([SAMTOOLS_EXE + ' faidx ' + filename + " " + query], stdout=PIPE, shell=True)
	seq = ""
	lines = process.stdout.readlines()
	if len(lines) > 1:
		for line in lines[1:]:
			seq += line.decode('ascii').strip()
	return seq.upper()


def sort_dataframe_by_score(df):
	df.sort(columns="CRISTA score", axis=0, ascending=False, inplace=True, kind='quicksort', na_position='last')
	df.reset_index(inplace=True, drop=True)
	df.index += 1 #start with 1


def remove_duplicates_from_df(df):
	sort_dataframe_by_score(df)
	df.drop_duplicates(subset=[CHROMOSOME_HEADING, STRAND_HEADING, STARTPOS_HEADING], inplace=True)
	df.drop_duplicates(subset=[CHROMOSOME_HEADING, STRAND_HEADING, ENDPOS_HEADING], inplace=True)


def readlines(f, newline):
	"""
	:param f: file or filepath
	:param newline: delimiter
	:return: a generator that reads the file as if line by line according to the delimiter

	:usage:
		with open('file') as f:
			for line in readlines(f, "."):
				print line
	"""
	if isinstance(f, str):
		f = open(f, "r")
	buf = ""
	while True:
		while newline in buf:
			pos = buf.index(newline)
			yield buf[:pos]
			buf = buf[pos + len(newline):]
		chunk = f.read(4096)
		if not chunk:
			yield buf
			break
		buf += chunk


def reverse_complement(seq):
	return Seq.reverse_complement(seq).upper()


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def get_bwa_genome_path(genome_assembly):
	dir = GENOMES_PATH + genome_assembly
	genome_path = glob.glob(dir + "/*.fa*.gz")[0]
	return genome_path


def parse_sam_out_table(samtable_filepath, sgrna_seq, genomic_db="hg19"):

	sam_gen = readlines(samtable_filepath, ";")

	first_line = sam_gen.__next__()
	firsts = first_line.split("XA:Z:")
	(first_line, alter) = firsts[0], None if len(firsts)==1 else firsts[1]

	data_df = pd.DataFrame(columns=[CHROMOSOME_HEADING, STRAND_HEADING, STARTPOS_HEADING,
									ENDPOS_HEADING, SGRNA_HEADING, OFFTARGET_HEADING])
	#parse best match (first line)
	first_line_tokens = first_line.split("\t")
	chromosome_num = "chr" + first_line_tokens[2]
	start_pos = int(first_line_tokens[3])
	cigar_str = first_line_tokens[5]
	dna_seq = first_line_tokens[9]
	
	if chromosome_num == "chr*" and start_pos == 0:
		return None, None

	features_mat = None
	is_first_line = True

	while alter and not alter.isspace() or is_first_line:
		#parse token
		if not is_first_line:
			alter_tokens = alter.split(",")
			chromosome_num = alter_tokens[0]

			chromosome_num = "chr" + chromosome_num
			strand = alter_tokens[1][0]
			start_pos = int(alter_tokens[1][1:])
			cigar_str = alter_tokens[2]

		# fetch sequence
		cigar_tokens = re.findall("([0-9]+[MIDNSHPX\=])", cigar_str)
		ot_pattern = ""
		for token in cigar_tokens:
			ot_pattern += (token[-1]*int(token[:-1]))
			ot_len = ot_pattern.count("M") + ot_pattern.count("S") + ot_pattern.count("D")

		end_pos = start_pos + ot_len - 1

		extended100_genomic_seq = fetch_dna_coordinates_offline(genomic_db, chromosome_num,
																start_pos - EXTENSION,
																end_pos + EXTENSION)
		#if "N" in extended100_genomic_seq:	# last change in 29/01/2018
		if "N" in extended100_genomic_seq or len(extended100_genomic_seq) < end_pos - start_pos + 1 + 2 * EXTENSION:
			if not is_first_line:
				alter = sam_gen.__next__()
			is_first_line = False
			continue

		if is_first_line:
			fwd_fuzzy_match = re.search("(" + dna_seq + "){i<=3,d<=3,s<=5,e<8}", sgrna_seq, re.I|re.BESTMATCH)
			rev_fuzzy_match = re.search("(" + reverse_complement(dna_seq) + "){i<=3,d<=3,s<=5,e<8}", sgrna_seq, re.I|re.BESTMATCH)
			fwd_fuzzy_cnt = float("inf") if fwd_fuzzy_match is None else sum(fwd_fuzzy_match.fuzzy_counts)
			rev_fuzzy_cnt = float("inf") if rev_fuzzy_match is None else sum(rev_fuzzy_match.fuzzy_counts)

			strand = "+" if fwd_fuzzy_cnt < rev_fuzzy_cnt else "-"

		extended100_genomic_seq = get_seq_by_orientation(extended100_genomic_seq, strand)
		dna_seq = extended100_genomic_seq[EXTENSION:-1 * EXTENSION]

		PAM_chunk = dna_seq[-2:]
		if PAM_chunk.count("G") + PAM_chunk.count("A") < 2 or "random" in chromosome_num:
			if not is_first_line:
				alter = sam_gen.__next__()
			is_first_line = False
			continue

		alignmentA, alignmentB, score = PA_script.align_pair(seqA=sgrna_seq[:-3], seqB=dna_seq[:-3], match_score=MATCH_SCORE,
															 mismatch_score=MISMATCH_PENALTY, gap_score=GAP_PENALTY,
															 gaps_allowed=3)
		alignmentA += sgrna_seq[-3:]
		alignmentB += dna_seq[-3:]

		features = CRISTA_main.get_features(aligned_sgRNA=alignmentA, aligned_offtarget=alignmentB, full_dna_seq=extended100_genomic_seq, pa_score=score)

		if features_mat is None:
			features_mat = np.asmatrix(features)
		else:
			features_mat = np.concatenate((features_mat, np.asmatrix(features)))

		data_df.loc[data_df.shape[0]] = [chromosome_num[3:], strand, start_pos, end_pos, alignmentA, alignmentB]

		if not is_first_line:
			alter = sam_gen.__next__()

		is_first_line = False

	return data_df, features_mat


def run_bwa(sgRNA_seq, genomic_db="hg19", n_gaps=N_GAPS,
			mm_penalty=MISMATCH_PENALTY*(-1), gap_open_penalty=GAP_PENALTY*(-1), gap_extension_penalty=0):
	"""
	:param sgRNA_seq: 23nt
	:param ref_genome_path:
	:param n_mms:
	:param n_gaps:
	:return: bwa output table for
	"""
	ref_genome_path = get_bwa_genome_path(genomic_db)
	sgrna_file = RESULTS_PATH + "sgrna.fa"
	sai_output_file = sgrna_file + ".out.sai"
	sam_output_file = sgrna_file + ".out.sam"
	output_table_file = sgrna_file + ".out.table"

	with open(sgrna_file, "w") as fpw:
		fpw.write(">seq\n" + sgRNA_seq)
	cmd = BWA_COMMAND_ON_SERVER + " aln "
	cmd += "-N " #don't stop once reaches the optimal match
	cmd += "-l 20 " #seed length
	cmd += "-i 0 " #Disallow an indel within INT bp towards the ends
	#cmd += "-e 3 " #Maximum number of gap extensions
	cmd += "-n {0} ".format(str(N_DISTEDIT+1)) #max number of mms in the whole region
	cmd += "-o {0} ".format(str(n_gaps)) #max number of gaps in the whole region
	cmd += "-d 3 " #Disallow a long deletion within INT bp towards the 3’-end
	cmd += "-k {0} ".format(str(N_DISTEDIT)) #Maximum edit distance in the seed
	cmd += "-M {0} ".format(str(0))
	cmd += "-O {0} ".format(str(1000000))
	cmd += "-E {0} ".format(str(gap_extension_penalty))
	cmd += ref_genome_path + " " + sgrna_file + " > " + sai_output_file
	os.system(cmd)
	if not os.path.exists(sai_output_file):
		raise Exception(IOError, "BWA samse failed, sai file not created")
	time.sleep(1)
	#logger.info("BWA align - execution of bwa completed")

	# convert sai to sam
	os.system(BWA_COMMAND_ON_SERVER + " samse -n 1000000 {0} {1} {2} > {3}".format(ref_genome_path, sai_output_file, sgrna_file, sam_output_file))
	if not os.path.exists(sam_output_file):
		raise Exception(IOError, "BWA samse failed, sam file not created")
	time.sleep(1)
	#logger.info("BWA samse - conversion from sai to sam completed!")

	# read sam file with samtools
	os.system(SAMTOOLS_EXE + " view -S " + sam_output_file + " > " + output_table_file)
	if not os.path.exists(output_table_file):
		raise Exception(IOError, "BWA samse failed, output table file not created")
	time.sleep(1)
	#logger.info(SAMTOOLS_EXE + " view - readable sam file conversion completed")

	data_df, features_mat = parse_sam_out_table(output_table_file, sgRNA_seq, genomic_db=genomic_db)
	#if data_df is None: # last change in 29/01/2018
	if features_mat is None:
		return None
	#data_df.to_csv(RESULTS_PATH + "data.csv")

	#score_df = CRISTA_main.predict_crista_score(features_mat)
	#logger.info("CRISTA predictions completed")
	#results_df = pd.concat((data_df, score_df), axis=1)

	#remove_duplicates_from_df(results_df)
	return data_df


def main(sgrna_seq, output_path, genome_assembly="TAIR10"):
	global RESULTS_PATH
	RESULTS_PATH = output_path

	### validate input: sgrna, genomic
	try:
		sgRNA_seq_re = re.search("[acgtu]+", sgrna_seq, re.IGNORECASE)
		assert sgRNA_seq_re is not None and len(sgRNA_seq_re.group()) == 20
	except:
		print("Invalid arguments. sgRNA sequence must be 20-nt long sequence of ACGTU nucleotides")
		print(parser.parse_args(['-h']))
		exit()

	sgrna_seq = sgrna_seq.upper() + "NGG"

	assert(os.path.exists(output_path), "Output path does not exist.")

	try:
		results_df = run_bwa(sgrna_seq, genomic_db=genome_assembly)
		if results_df is None:
			raise ValueError("Matches were not found in the selected genome reference.")
		filename = sgrna_seq + "_CRISTA_offtargets_scores.csv"
		results_df.to_csv(RESULTS_PATH + "/" + filename)

	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)


if __name__ == '__main__':

	logger = logging.getLogger('CRISTA main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='Looks for off-targets for a single sgRNA')
	parser.add_argument('--sgseq', '-s', default=None, help="20bp long (without PAM)", required=True)
	parser.add_argument('--genome_database', '-g', default="arabidopsis_thaliana")
	parser.add_argument('--output_path', '-p', default=None, required=True)

	args = parser.parse_args()

	sgrna_seq = args.sgseq
	genome_assembly = args.genome_database
	output_path = args.output_path

	main(sgrna_seq, output_path, genome_assembly)
