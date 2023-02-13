import os
import re
import pickle
import numpy as np
from cfd_score_calculator import *


# Param1: file_path- CRISTA output file
# Param2: optional- accounting for a specific sgRNA or not. The default is None (means not accounting for)
# Return value: dictionary CRISTA data and cfd score
def parse_crista_results(file_path, sgRNA = None):
    dic = {}
    file = open(file_path, 'r')
    file.readline()
    for line in file:
        gene_data = re.split(",", line)
        p = re.compile("[\S]+")
        digit_pattern = re.compile("[\d]+")
        gene_index =  p.findall(gene_data[0])[0]
        chromosome = p.findall(gene_data[1])[0]
        strand =  p.findall(gene_data[2])[0]
        start_pos = p.findall(gene_data[3])[0]
        end_pos = p.findall(gene_data[4])[0]

        sgrnas = get_input_for_cfd(p.findall(gene_data[5])[0], 1)#, chromosome, strand, int(float(start_pos)), int(float(end_pos)))
        site =  get_input_for_cfd(p.findall(gene_data[6])[0], 0)[0]#, chromosome, strand, int(float(start_pos)), int(float(end_pos)))[0]
        score = str(calculate_cfd_score(sgrnas, site))
        if score == '-1':
            continue    #If the alignment is gapped, It should be ignored, and not put to the dictionary
        if digit_pattern.findall(chromosome) != []: #if not a digit, then it is mitochondria or chloroplast
            if not sgRNA:
                dic[gene_index] = (int(chromosome),strand, int(float(start_pos)), int(float(end_pos)), float(score))
            else:
                dic[gene_index] = (int(chromosome), strand, int(float(start_pos)), int(float(end_pos)), float(score)
                                   , sgRNA, site[:20])
        elif chromosome == "chloroplast" or chromosome == 'Pt': # set chloroplast as chromosome -1
            if not sgRNA:
                dic[gene_index] = (-1,strand, int(float(start_pos)), int(float(end_pos)), float(score))
            else:
                dic[gene_index] = (-1, strand, int(float(start_pos)), int(float(end_pos)),
                                   float(score), sgRNA, site[:20])

        else:   # it is mitochondria- set mitochondria as -2
            if not sgRNA:
                dic[gene_index] = (-2, strand, int(float(start_pos)), int(float(end_pos)), float(score))
            else:
                dic[gene_index] = (-2, strand, int(float(start_pos)), int(float(end_pos)), float(score),
                                   sgRNA, site[:20])
    file.close()
    return dic



#######################################################################################################

# Param1: sgrnas- a list possible variations of sgRNAs (replacing N in the pam segment by {A, T, C, G}
# Param2: aligned site of the sgRNA.
# Return value: cfd score


def calculate_cfd_score(sgrnas, site):
    #print("site is : "+site+"\n")
    #print("sgRNAs are: "+str(sgrnas)+"\n")

    dic_nuc = {'N': {'T', 'A', 'G', 'C'}, 'Y':{'T', 'C'}, 'R':{'A', 'G'}, 'K':{'G', 'T'}, 'M':{'A', 'C'},
               'S':{'G', 'C'}, 'W':{'A', 'T'}, 'D':{'T', 'A', 'G'}}
    max_score = 0

    for sgrna in sgrnas:
        new_aligned_site = ""
        for i in range(len(site)):
            if site[i] in dic_nuc:
                char_in_sg = sgrna[i]
                for char in dic_nuc[site[i]]:
                    if char == char_in_sg:
                        new_aligned_site += char

                if len(new_aligned_site) == i:
                    new_aligned_site += list(dic_nuc[site[i]])[0]

            else:
                new_aligned_site += site[i]


        #print(sgrna+"\n")
        #print(new_aligned_site+"\n")
        score = calculate_cfd(sgrna[:20], new_aligned_site[:20])
        if score == None:
            #print("sgrna is: "+sgrna+"\n")
            #print("aligned region is: "+ new_aligned_site+"\n")
            return -1 #TODO when does it happen????
        if score > max_score:
            max_score = score
    #if sgrna[:20] == "ATGTCCACGTCTTAAAGTTT" and new_aligned_site[:20] == "CTGTTCACATCTTAATGTTT":
        #print("sgRNA is : "+ sgrna[:20] + " target is : "+ new_aligned_site[:20] +" and score is "+ str(max_score)+"\n")
    return max_score

######################################################################################################################


def get_input_for_cfd(seq, sgrna_bool):
    new_sequences = []
    nucleotides = ['A', 'T', 'C', 'G']
    if sgrna_bool:
        for nucleotide in nucleotides:
            new_sequences.append(re.sub('N', nucleotide, seq))
    else:
        new_sequences.append(seq)
    return new_sequences


###################################################################################################
#input: annotation file path, destination pickle file
#output: pickle file
# No return value

def parse_annotation_file(source_path, pickle_path):

    dic = {}
    file = open(source_path, 'r')
    file.readline()
    p = re.compile("[\S]+")

    for line in file:
        gene_data = extract_data_annotation(line)
        splitted_line = re.split(";", line)
        gene_id = p.findall(splitted_line[0])[0][1:-1]

        dic[gene_id] = gene_data

    file.close()
    pickle.dump(dic, open(pickle_path, "wb"))
    return



###############################################################################
#A function which finds the file path of the crista output file given the parent directory
#Return: full path of the CRISTA output file.


def find_crista_output_filename(search_directory, sgRNA):
    file_name = os.path.join(search_directory, sgRNA +"NGG"+ "_CRISTA_offtargets_scores.csv")
    if not os.path.exists(file_name):
        return None
    return file_name
####################################################################################################
# Param1: gene_whole_data- a line of an entry in annotation file
# Return value: gene's data, which includes chromosome, strand ['+'|'-'], exon positions, start cds, stop cds,
# start transcript, stop transcript

def extract_data_annotation(gene_whole_data):

    p = re.compile("[\S]+")
    pattern = re.compile("[\d]+")

    splitted_line = re.split(";", gene_whole_data)
    start = pattern.findall(splitted_line[4])[0]
    stop = pattern.findall(splitted_line[5])[0]
    start_transcript, stop_transcript = find_transcript_ends(p.findall(splitted_line[6])[0])
    exon_positions = extract_exons_positions(p.findall(splitted_line[3])[0])
    strand = p.findall(splitted_line[7])[0][1:-1]
    chromosome_raw = pattern.findall(splitted_line[8])

    chromosome = determine_chromosome(chromosome_raw, splitted_line)

    return (int(chromosome), strand, exon_positions, int(start), int(stop), int(start_transcript), int(stop_transcript))

####################################################################################################################33


def find_transcript_ends(transcript_fragments):
    p = re.compile("([\d]+)..([\d]+)")
    transcript_fragments_str = p.findall(transcript_fragments)
    start = int(transcript_fragments_str[0][0])
    end = int(transcript_fragments_str[-1][1])
    return start, end

#######################################################################################################################


# Param1: chromosome_raw-> string representation of chromosome
# Param2: splitted line-> a splitted list of all the gene's data
# Return value: chromosome in digital representation


def determine_chromosome(chromosome_raw, splitted_line):
    p = re.compile("[\S]+")

    chromosome = None

    if chromosome_raw != []:
        chromosome = chromosome_raw[0]
    else:
        chromosome_raw = p.findall(splitted_line[8])[0][1:-1]
        if chromosome_raw == "chloroplast":
            chromosome = -1
        elif chromosome_raw == "mitochondrial":
            chromosome = -2
    return chromosome
#######################################################################################################

# This function creates a dictionary, where the key is the chromosome number and the value is the its length
# Param1: file path of the chromosomes' sequences
# Return value: dictionary where the key is the chromosome number and the value is its sequence.
# Note: Mitochondrial genome is denoted by -2, and the chloroplast genome is indicated by -1.

def create_chromosome_seq_file(file_path):
    source_file = open(file_path, 'r')
    content = source_file.read()
    source_file.close()
    delimiter_chr = re.compile(">([\S]+)([^>]+)")


    seqs_chromosomes = delimiter_chr.findall(content)
    dic = {}
    for i in range(len(seqs_chromosomes)):
        splitted_chr_seq = re.split("\n|[\s]+", seqs_chromosomes[i][1])
        sequence = "".join(splitted_chr_seq)
        #dic[i+1] = len(sequence)
        if seqs_chromosomes[i][0] == 'chloroplast':
            dic[-1] = len(sequence)
        elif seqs_chromosomes[i][0] == 'mitochondrial':
            dic[-2] = len(sequence)
        else:
            dic[int(seqs_chromosomes[i][0])] = len(sequence)

    return dic

########################################################################################################

# Param1: exon postions in a string
# Return value: exon positions list, where each exon is represented as a tuple (start index, stop index)

def extract_exons_positions(exons_positions_str):
    p = re.compile("([\d]+)..([\d]+)")
    positions_tuples_str = p.findall(exons_positions_str)
    final_exon_postions = []
    for position_range in positions_tuples_str:
        final_exon_postions.append((int(position_range[0]), int(position_range[1])))
    return final_exon_postions

######################################################################################3
# return a dictionary where the key is chromosome, and the value is a list of [gene exons, strand, gene]
def getExonsDictionary(annotation_data):
    dict_exons = {}
    for gene in annotation_data:
        chromosome = annotation_data[gene][0]
        strand = annotation_data[gene][1]
        exons = annotation_data[gene][2]
        if chromosome in dict_exons:
            dict_exons[chromosome].append([exons, strand, gene])
        else:
            dict_exons[chromosome] = [[exons, strand, gene]]
    return dict_exons


