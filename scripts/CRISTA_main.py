#########################################################################
#########################################################################
##                                                                     ##
##                                                                     ##
##                              CRISTA                                 ##
##                                                                     ##
##                                                                     ##
##                 A tool for CRISPR Targets Assessment                ##
##                                                                     ##
##                             v. 1.0                                  ##
##                                                                     ##
##                                                                     ##
#########################################################################
#########################################################################
## This code is provided by Shiran Abadi                               ##
##                                                                     ##
## CRISTA is based on learning a regression model using the Random     ##
## Forest algorithm within the machine learning paradigm. CRISTA can   ##
## be used to determine the propensity of a genomic site to be cleaved ##
## by a given sgRNA. CRISTA was trained on a large dataset assembled   ##
## from published data of genome-wide unbiased methods for CRISPR-Cas9 ##
## cleavage sites profiling [1–5]. It accounts for the possibility of  ##
## bulges and incorporates a wide range of features encompassing those ##
## that are specific to the genomic content, features that define the  ##
## thermodynamics of the sgRNA, and features concerning the pairwise   ##
## similarity between the sgRNA and the genomic target. Altogether,    ##
## these form a complex model that can be used to predict the          ##
## cleavage propensity of a selected genomic site.                     ##
##                                                                     ##
## More functionalities are available at www.crista.tau.ac.il          ##
##                                                                     ##
## For academic use, please cite crista.tau.ac.il.                     ##
## Non-commercial use!                                                 ##
##                                                                     ##
## Please do not change and distribute.                                ##
##                                                                     ##
#########################################################################
#########################################################################
##                                                                     ##
##    usage: command line                                              ##
##    python CRISTA.py -s SGRNA_SEQ -d GENOMIC_SEQ                     ##
##                                                                     ##
##                                                                     ##
##     SGRNA_SEQ: sgRNA sequence of 20 bases (without PAM)             ##
##     GENOMIC_SEQ: DNA target sequence with 3 additional bases at     ##
##                   each end (total of 29 nucleotides)                ##
##                                                                     ##
#########################################################################
#########################################################################
##                                                                     ##
##    Dependencies:                                                    ##
##       python 3                                                      ##
##       numpy, sklearn, pickle, and argparse modules                  ##
##                                                                     ##
#########################################################################
# usage example
# python CRISTA.py -s CTCAGCTGAGGTTGCTGCTG -d GGCCTCAGCTGAGGTTGCTGCTGTGGAAG
#########################################################################

import argparse
import pickle
import random
import re
import numpy as np
import pandas as pd

### globals
RF_MODEL_DIR_WOGENOMIC = "/bioseq/crista/rf_nogenomic/"
MATCH_SCORE = 1.0
MISMATCH_PENALTY = 0.0
GAP_PENALTY = -1.25
MAX_ALLOWED_GAPS = 3
EXTENSION = 3
DNA_PAIRS_THERMODYNAMICS = {"AA": 9.1, "AT": 8.6, "TA": 6.0, "CA": 5.8, "GT": 6.5, "CT": 7.8, "GA": 5.6, "CG": 11.9,
							"GC": 11.1, "GG": 11.0, "TT": 9.1, "TG": 5.8, "AC": 6.5, "AG": 7.8, "TC": 5.6, "CC": 11.0} #Breslauer et al.
DNASHAPE_DICT_FILE = "/groups/itay_mayrose/shiranabad/CRISPR/data/dnaShape.pkl"
DNASHAPE_DICT = None
ACGT_REPLACEMENT = {"A": '1', "C": '2', "G": '3', "T": '4', 'N': '0'}

MMS_TYPE_REPLACEMENT = {'0': "match", '-1': "indel", '1': "wobble", '2': "RR transition", '3': "YY transition", '4': "transversion"}


def agct2numerals(st):
	new_st = ""
	for x in st:
		try:
			new_st += ACGT_REPLACEMENT[x]
		except KeyError:
			new_st += "0"
	return new_st


def get_avg(l):
	return sum(l)/float(len(l))


def count_mismatches(aligned_seq1, aligned_seq2):
	"""
	:param seq1, aligned_seq2: aligned sgRNA and genomic target (seq+PAM)
	:return:
	"""
	cnt = 0
	ending = len(aligned_seq1) - 3

	for i in range(ending):
		cnt += int(aligned_seq1[i] != aligned_seq2[i] and aligned_seq1[i] != "-" and aligned_seq2[i] != "-")
	return cnt


def cnt_bulge(aligned_seq):
	return aligned_seq.count("-")


def count_consecutive_inconsistencies(aligned_seq1, aligned_seq2):
	"""
	:param seq1, aligned_seq2: aligned sgRNA and genomic target (seq+PAM)
	:return: number of concatenated-extended mismatches and bulges
	"""
	cnt = 0
	current_cnt = 0

	for i in range(len(aligned_seq2) - 3):
		if aligned_seq2[i] != aligned_seq1[i]:
			current_cnt += 1
		else:
			cnt += current_cnt > 0
			current_cnt = 0
	return cnt


def get_DNA_enthalpy(dna_seq):
	seq = re.sub("[^ACGT]", lambda x: random.choice(["A", "C", "G", "T"]), dna_seq)
	return sum([DNA_PAIRS_THERMODYNAMICS[seq[i-1:i+1]] for i in range(1, len(seq))])


def get_DNAshape_features(dna_seq):
	"""
	:param dna_seq: sequence of nucleotides
	:return: a dictionary with scores of rigidity for Major Groove Width (MGW), ProT (Propeller-Twist), Roll, and HelT (Helical-Twist).
	The values are the scores for each pentamer/hexamer as computed by DNAshape (Zhou et al., doi:10.1093/nar/gkt437)
	across the DNA sequence
	"""

	global DNASHAPE_DICT
	if DNASHAPE_DICT is None:
		DNASHAPE_DICT = pickle.load(open(DNASHAPE_DICT_FILE, "rb"))

	mgw = [None]
	roll = [None]
	prot = [None]
	helt = [None]

	for i in range(2, len(dna_seq)-2):
		current_heptamer = dna_seq[i-2 : i+3]
		current_heptamer = re.sub("[^ACGT]", lambda x: random.choice(["A", "C", "G", "T"]), current_heptamer)
		current_nucleotide = DNASHAPE_DICT[current_heptamer]
		mgw += current_nucleotide["MGW"]
		roll += current_nucleotide["Roll"]
		prot += current_nucleotide["ProT"]
		helt += current_nucleotide["HelT"]

	helt_modified = [helt[1]]
	for i in range(2, len(helt), 2):
		helt_modified.append(get_avg(helt[i:i+2]))
	roll_modified = [roll[1]]
	for i in range(2, len(roll), 2):
		roll_modified.append(get_avg(roll[i:i+2]))

	return {"MGW":mgw[1:], "ProT": prot[1:], "Roll": roll_modified, "HelT": helt_modified}


def get_features(full_dna_seq, aligned_sgRNA, aligned_offtarget, pa_score):
	"""
	compute CRISTA features
	"""

	# get alignment features
	mms_cnt = count_mismatches(aligned_sgRNA, aligned_offtarget)
	rna_bulges = cnt_bulge(aligned_sgRNA)
	dna_bulges = cnt_bulge(aligned_offtarget)
	gapless_dnaseq = re.sub("-", "", aligned_offtarget)

	# quartets mismatches counts
	rev_rna = (aligned_sgRNA[::-1])[3:]
	rev_dna = (aligned_offtarget[::-1])[3:]
	mismatches_1_4 = count_mismatches(rev_rna[:4], rev_dna[:4])
	mismatches_5_8 = count_mismatches(rev_rna[4:8], rev_dna[4:8])
	mismatches_9_12 = count_mismatches(rev_rna[8:12], rev_dna[8:12])
	mismatches_13_16 = count_mismatches(rev_rna[12:16], rev_dna[12:16])
	mismatches_17_end = count_mismatches(rev_rna[16:], rev_dna[16:])

	# get from alignment mismatches per position
	mismatches = [-2] * 23  # 5' -> 3', without PAM                                             # undefined
	offset = 26 - len(aligned_offtarget)
	for i in range(len(aligned_offtarget) - 3):
		rna_base = aligned_sgRNA[i]
		dna_base = aligned_offtarget[i]
		# Categorization of mismatch type
		if rna_base == dna_base:                                                                # match
			mismatches[offset + i] = 0
		elif rna_base == "-" or dna_base == "-":                                                # indel
			mismatches[offset + i] = -1
		elif (rna_base == "T" and dna_base == "C") or (rna_base == "G" and dna_base == "A"):    # wobble: rG:dA, rT:dC
			mismatches[offset + i] = 1
		elif rna_base in ["G", "A"] and dna_base in ["T", "C"]:                                 # R-R pairing
			mismatches[offset + i] = 2
		elif rna_base in ["C", "T"] and dna_base in ["G", "A"]:                                 # Y-Y pairing
			mismatches[offset + i] = 3
		elif (rna_base == "A" and dna_base == "G") or (rna_base == "C" and dna_base == "T"):    # other transversion
			mismatches[offset + i] = 4

	# total types mismatches
	wobble_total = mismatches.count(1)
	RR_total = mismatches.count(2)
	YY_total = mismatches.count(3)
	Tv_total = mismatches.count(4)

	# pairs of nucleotides in positions 1-5 upstream to PAM (1-2, 2-3, 3-4, 4-5)
	seed_couples = []
	for i in range(4):
		seed_couples.append(agct2numerals(gapless_dnaseq[-8 + i:-6 + i]))

	# PAM and 5'-end nucleotides
	PAM_2_first = agct2numerals(gapless_dnaseq[-2:])
	PAM_N_id = agct2numerals(gapless_dnaseq[-3])
	last_pos_nucleotide = agct2numerals(gapless_dnaseq[0])

	# mismatches and bulges - linked
	consecutive_inconsistencies_cnt = count_consecutive_inconsistencies(aligned_sgRNA, aligned_offtarget)
	avg_inconsistency_length = (mms_cnt + rna_bulges + dna_bulges) / float(
		consecutive_inconsistencies_cnt) if consecutive_inconsistencies_cnt > 0 else 0

	# nucleotides occupancies in DNA target sequence
	nA = gapless_dnaseq.count("A")
	nC = gapless_dnaseq.count("C")
	nG = gapless_dnaseq.count("G")
	nT = gapless_dnaseq.count("T")

	# GC content
	extended_genomic_GC_content = (full_dna_seq.count("C") + full_dna_seq.count("G")) / float(
		len(full_dna_seq))

	# five nucleotides downstream to PAM
	dna_plus_3 = full_dna_seq[EXTENSION - 3 :EXTENSION + len(gapless_dnaseq) + 6]
	nucleotides_down_pam = [agct2numerals(c) for c in full_dna_seq[EXTENSION + len(gapless_dnaseq) + 1:EXTENSION + len(
		gapless_dnaseq) + 6]]  # the model feature for additional two nucleotides is disregarded (0) but still exists

	# geometry features: dna_enthalpy
	extended_dna_enthalpy = get_DNA_enthalpy(full_dna_seq)
	dna_enthalpy = get_DNA_enthalpy(gapless_dnaseq)

	# geometry features: DNA shape per pentamer
	dna_shape_features = get_DNAshape_features(dna_plus_3)

	features = [aligned_sgRNA, aligned_offtarget,
				pa_score, rna_bulges + dna_bulges, rna_bulges, dna_bulges, PAM_2_first,
	        PAM_N_id, last_pos_nucleotide, mms_cnt,
	        consecutive_inconsistencies_cnt, avg_inconsistency_length] \
	       + [mismatches_1_4, mismatches_5_8, mismatches_9_12, mismatches_13_16, mismatches_17_end] + \
	       [wobble_total, YY_total, RR_total, Tv_total] + \
	       seed_couples + [agct2numerals(gapless_dnaseq[i]) for i in range(-8, -3)] + \
	       [extended_genomic_GC_content, #upstream_50_extension_gc, downstream_50_extension_gc,
	        dna_enthalpy, extended_dna_enthalpy,
	        nA, nC, nT, nG] + nucleotides_down_pam + \
	       [min(dna_shape_features["MGW"])] + \
	       [get_avg(dna_shape_features["HelT"])] + \
	       [get_avg(dna_shape_features["Roll"])] + \
	       [get_avg(dna_shape_features["ProT"])] + \
	       [dna_shape_features["MGW"][-3], dna_shape_features["HelT"][-3],
	        dna_shape_features["Roll"][-3], dna_shape_features["ProT"][-3]]

	return np.array(features).reshape((1,len(features)))


def predict_crista_score(features_mat, logger=None):
	"""
	:param features_df: dataframe: first col: rna, second: dna, the rest are features
	mode: either full, nogenomic or noflanking
	:return: features df + prediction col
	"""
	n_predictors = 10 #100

	path = RF_MODEL_DIR_WOGENOMIC

	predictors = [path + str(i)+ "/RFpredictor.pkl" for i in range(n_predictors)]
	prediction_df = pd.DataFrame()

	features_mat = np.delete(features_mat, [0,1], axis=1)

	for i in range(n_predictors):
		rf_predictor = pickle.load(open(predictors[i], "rb"))
		prediction_df[i] = rf_predictor.predict(features_mat)
		if i % 10 == 9 and logger:
			logger.info(str(i + 1) + "% are done")
	final_pred_df = pd.DataFrame(np.divide(prediction_df.mean(axis=1), 8.22), columns=["CRISTA score"])
	final_pred_df.ix[final_pred_df["CRISTA score"]>1, "CRISTA score"] = 1.0

	return final_pred_df