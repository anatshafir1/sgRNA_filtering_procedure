# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:48:30 2017

@author: User
"""
#from create_dictionaries_and_pickle import *
import numpy as np
import re
import os
from global_variables import *
########################################################################################################################
# This function checks whether a given sgRNA is specific to on-target genes. If it has a particular high scored off-target
# in the exon zone, it should be filtered out.
# Param1: threshold in the range (0,1) which is a criterion for the sgRNA to have a high or low affinity to possible
# off-targets. Could be either relative to the max on-target score or not.
# Param2: True if the threshold is relative to the max on-target score. False- otherwise.
# Param3: a list of on-target scores (the scores for the genes found in Crispys).
# Param4: dictionary which contains the information of the targets found by BWA:
# key = target index, value = (chromosome,strand, start pos, stop pos, score)
# Param5: a path for the annotation file of the family
# Param6: dictionary of chromsomes' lengths
# Return Value: True if the sgRNAs should not be filtered. False if it should be filtered.
def check_sgRNAs_compatability_alternative(threshold, relative, on_target_scores, dict_of_CRISTA_results,
                                           annotation_file_path, dic_genome):
    path_for_exons_dic = os.path.join(MAIN_FILES_DIR, 'dict_exons.p')
    with open(path_for_exons_dic, 'rb') as handle:
        dict_exons = pickle.load(handle)
    handle.close()
    with open(annotation_file_path, 'rb') as handle:
        annotation_data = pickle.load(handle)
    handle.close()
    for gene_index in dict_of_CRISTA_results:
        score = dict_of_CRISTA_results[gene_index][-1]
        high_score = check_if_score_is_high(score, relative, on_target_scores, threshold)
        if not high_score:
            continue
        on_target = get_ontarget_gene(gene_index, dict_of_CRISTA_results, annotation_data, dic_genome)
        if on_target:
            continue
        # consider the option of off- target
        chromosome = dict_of_CRISTA_results[gene_index][0]
        # search only in the relevant chromosome
        exons_and_data = dict_exons[chromosome]

        for exon_data in exons_and_data:
            gene = exon_data[2]
            exons = exon_data[0]
            strand = exon_data[1]
            if gene in annotation_data: #on- target but the cut zone is not in exon
                continue
            # check if off-target (and if it is cleaved in the exon zone)
            off_target = check_if_off_target(gene_index, dict_of_CRISTA_results, exons, strand, dic_genome[chromosome])
            if off_target:
                return False
    return True

####################################################################################################################
# checks if a particular target is within an off-target gene (in one of its exons)
def check_if_off_target(gene_index, dic_crista, exons, strand, chromosome_length):
    off_target = False
    strand_crista = dic_crista[gene_index][1]
    stop = dic_crista[gene_index][3]
    nucleotides_of_cut = [stop-3-5+1, stop-3-4+1]   # 5 nucleotides upstream or 4 nucleotides upstream PAM
    if strand_crista != strand:
        nucleotide_of_cut_in_gene1 = chromosome_length - nucleotides_of_cut[0] + 1
        nucleotide_of_cut_in_gene2 = chromosome_length - nucleotides_of_cut[1] + 1
    else:
        nucleotide_of_cut_in_gene1 = nucleotides_of_cut[0]
        nucleotide_of_cut_in_gene2 = nucleotides_of_cut[1]
    for exon in exons:
        start_exon = exon[0]
        stop_exon = exon[1]
        if (start_exon <= nucleotide_of_cut_in_gene1 <= stop_exon) or\
                (start_exon <= nucleotide_of_cut_in_gene2 <= stop_exon):
            off_target = True
            break
    return off_target


########################################################################################################################
# check if a given off-target has a excessively high score
def check_if_score_is_high(score, relative, on_target_scores, threshold):
    if relative:
        return score > float(threshold)*max(on_target_scores)

    return score > float(threshold)

########################################################################################################################

## Input ##
#Param1: gene id as it appears in CRISTA
# param 2: crista dictionary
# param3: all related family genes
# param4: annotation file dictionary, which
#contains gene name as a key, and all the related data as a value
# param5: dictionary of the whole genome of arabidopsis Thaliana-
#the key is the chromosome and the value is the its length; param5: a file path, the file will contain the binding genes,
#which are on target, and their location (intron or exon).
##Output##
#Return : gene name if the binding site of the sgRNA is ontarget (or None if it isn't) and a boolean value which will
#indicate if the binding site is located in intron or exon. If the binding site is off- target, the value will be False.
#It will be true only if the bs inside an exon.

def get_ontarget_gene(gene_index, dic_crista, annotation_data, dic_genome):

    family_gene_target = None
    chromosome = dic_crista[gene_index][0]
    strand = dic_crista[gene_index][1]
    #start = dic_crista[gene_index][2]
    stop = dic_crista[gene_index][3]
    nucleotides_of_cut = [stop - 3 - 5 + 1, stop - 3 - 4 + 1]  # 5 nucleotides upstream or 4 nucleotides upstream PAM

    for gene in annotation_data:  # checking if the gene is related to the family

        data_of_gene = annotation_data[gene]
        chromosome_annotation = data_of_gene[0]
        strand_annotation = data_of_gene[1]
        start_cds = data_of_gene[3]
        stop_cds = data_of_gene[4]


        if chromosome_annotation != chromosome:
            continue

        elif strand != strand_annotation:
            chromosome_length = dic_genome[chromosome]
            #chromosome_length = dic_genome[chromosome]
            nucleotide_of_cut_in_gene1 = chromosome_length - nucleotides_of_cut[0] + 1
            nucleotide_of_cut_in_gene2 = chromosome_length - nucleotides_of_cut[1] + 1

        else:
            nucleotide_of_cut_in_gene1 = nucleotides_of_cut[0]
            nucleotide_of_cut_in_gene2 = nucleotides_of_cut[1]

        if (stop_cds >= nucleotide_of_cut_in_gene1 >= start_cds ) or (stop_cds >= nucleotide_of_cut_in_gene2 >= start_cds):
            family_gene_target = gene
            return family_gene_target



    return family_gene_target




#######################################################################################################################
# Param1: A cluster of sgRNAs represented by a subset of genes (key = sorted subgroup of genes, value = list of sgRNAs)
# Param2: dictionary of cds of the family genes (key = gene, value = list of exons)
# Param3: list of sgRNAs the so far remained sgRNAs candidates.
# Pram4: a dictionary- key = sgRNA sequence, value = (set of genes, dictionary of positions, overall score)
# #   where dictionary of positions is: key = gene, value = list of positions where each element is a tuple (position, boolean),
# #   The boolean is True if the strand is + and Flase if the strand is -.
# Param5: threshold distance for overlap (a fraction of the positions in sgRNA)
# Param6: Crispys threshold (omega)
# Return Value: updated set of clustered sgRNAs per subgroup

def filter_same_positioned_sg(sets_of_sgRNAs, dict_of_cds, proccessed_sg_list, dict_of_pos,
                              thre_distance, crispys_threshold):
    copy_set_of_sgRNAs = {}
    removed_sgs = set()
    threshold_distance = float(thre_distance)
    dict_of_sg_candidate = create_dict_sg_candidate(proccessed_sg_list) # dict: key = sgRNA seq, value = candidate inst.
    dict_of_subsets, sorted_sets = create_dict_of_subsets(sets_of_sgRNAs)    #for each subset of genes, there is a list of its subsets

    for gene_set in sorted_sets:    #gene set is now the main subset
        subsets = dict_of_subsets[gene_set] # sorted (in a decending order) subsets of genes of the current subset(represented as sets)
        set_sgs_main_subset = sets_of_sgRNAs[gene_set]
        for i in range(len(subsets)):   #iterating over subsets
            subset_sorted = tuple(sorted(list(subsets[i])))
            set_sgs_subset = sets_of_sgRNAs[subset_sorted]   #set of sgRNAs
            for sgrna1 in set_sgs_main_subset:
                if sgrna1 in removed_sgs:
                    continue
                for sgrna2 in set_sgs_subset:
                    if sgrna2 in removed_sgs:
                        continue
                    if sgrna1 != sgrna2:
                        decide_whether_to_filter(dict_of_sg_candidate, crispys_threshold, dict_of_pos, dict_of_cds,
                                                 threshold_distance, removed_sgs, sgrna1, sgrna2, proccessed_sg_list)
    create_copy_set_of_sgs(sets_of_sgRNAs, removed_sgs, copy_set_of_sgRNAs)
    return copy_set_of_sgRNAs
#######################################################################################################################
# set of internal nodes represented as tuples of genes in a sorted order.

def create_set_of_internal_nodes(lst_sgRNAs):
    lst_of_internal_nodes = [tuple(sorted(subgroup.genes_lst)) for subgroup in lst_sgRNAs]
    set_of_internal_nodes = set(lst_of_internal_nodes)
    return set_of_internal_nodes
########################################################################################################################

def create_dict_sg_candidate(proccessed_sg_list):
    dic = {}
    for candidate in proccessed_sg_list:
        sg = candidate.seq
        dic[sg] = candidate
    return dic


########################################################################################################################
# Param1: file path of exons for each gene of the family.
# Return: a dictionary where the key is the gene, and the value is a list of exons.


def extract_cds(file_path):

    dic = {}
    p = re.compile(">([\S]+)[\s]+([\S]+)", re.MULTILINE)
    file = open(file_path, 'r')
    content = file.read()
    file.close()
    genename_and_seq = p.findall(content)
    for entry in genename_and_seq:
        if entry[0] not in dic:
            dic[entry[0]] = [entry[1]]
        else:
            dic[entry[0]].append(entry[1])
    return dic

##################################################################################################################

def create_dict_score_sgs(candidates):
    dic_scores = {}
    for candidate in candidates:
        dic_scores[candidate.seq] = sum(list(candidate.genes_score_dict.values()))
    return dic_scores

################################################################################################################

# dictionary with sgRNA seq as a key and the total score of genes above threshold as a value
def create_dict_score_threshold(proccessed_sg_list, threshold):
    dic_sg_threshold = {}
    for candidate in proccessed_sg_list:
        score_of_sg = 0
        dict_gene_score = candidate.genes_score_dict
        for gene in dict_gene_score:
            if dict_gene_score[gene] >= threshold:
                score_of_sg += dict_gene_score[gene]
        dic_sg_threshold[candidate.seq] = score_of_sg
    return dic_sg_threshold



#######################################################################################################################
# Param1: dictionary: key = subgroup, value = set of sgRNAs
# Param2: if candidate=1, the instance in subgroups_and_sgs is candidate, otherwise- 0.
# Param3: the required number of sgRNA per internal node
# Param4: a dictionary- key = sgRNA seq, value = the total score
# Param5: a set of subgroups of genes (as they appear in Crispys output).
# Param6: boolean- true if bottom up
# Param7: # update a dictionary dict_sgrna_score_threshold in the following manner: key = sgRNA seq, value= overall score of
# # sgRNA, composed only of genes above threshold
# Return value: dictionary, where key=subgroup of genes and value = best sgRNAs limited to the permitted number of sgRNAs.

def create_dict_of_sgs_per_subgroup_alternative(subgroups_and_sgs, candidate, num_of_sgs, dict_scores,
                                                internal_nodes, bottom_up_indicator, dict_sgrna_score_threshold):

    all_sets = set(subgroups_and_sgs.keys()).union(internal_nodes)
    dict_of_key_subsets = create_dict_with_dummy_values(all_sets)
    # sorting the subgroups in ascending order
    dict_of_subsets, sorted_sets = create_dict_of_subsets(dict_of_key_subsets, bottom_up=bottom_up_indicator,
                                                          top_down=not bottom_up_indicator)
    dict_sgs_per_subgroup = {}
    all_sgrnas = set() #this set will contain all the sgRNAs from all sets.

    for gene_set in sorted_sets:
        if not gene_set in subgroups_and_sgs:
            continue
        all_sgrnas = all_sgrnas.union(add_seq_to_gene_set(subgroups_and_sgs, gene_set, candidate))
    for gene_set in sorted_sets:
        if not gene_set in internal_nodes:
            continue
        list_of_sgs_to_add = find_potential_sg_to_add(all_sgrnas, gene_set, dict_scores, dict_of_subsets,
                                                      dict_sgrna_score_threshold,
                                                      alternative = (subgroups_and_sgs, candidate))
        dict_sgs_per_subgroup[gene_set] = set()
        for i in range(min(num_of_sgs, len(list_of_sgs_to_add))):
            # add the best sgRNAs, and remove them from the set of potential sgRNAs.
            dict_sgs_per_subgroup[gene_set].add(list_of_sgs_to_add[i])
            all_sgrnas -= {list_of_sgs_to_add[i]}
    return dict_sgs_per_subgroup
##################################################################################################################################################################
# Param1: all the available sgRNAs for the current subgroup of genes
# Param2: the current subgroup of genes
# Param3: a dictionary- key = sgRNA seq, value = the total score
# Param4: a dictionary of subgroups of genes and their subsets
# Param5: # update a dictionary dict_sgrna_score_threshold in the following manner: key = sgRNA seq, value= overall score of
# # # sgRNA, composed only of genes above threshold
# Param6: optional- alternative. If not none: (subgroups_and_sgs, candidate), where subgroups_and_sgs is a dictionary:
# key = subgroup, value = set of sgRNAs and candidate is 1 or 0 (candidate inst. or sgRNA seq).
# returns a sorted list of potential sgRNAs for a given subgroup of genes (sorted according to the score in a descending order).

def find_potential_sg_to_add(sgrna_for_next_round, gene_set, dict_scores, dict_of_subsets, dict_of_threshold_scores, alternative = None):
    set_of_sg_to_use = set()
    subsets = dict_of_subsets[gene_set]
    if not alternative:
        for subset in subsets:
            sorted_subset = tuple(sorted(list(subset)))
            if sorted_subset in sgrna_for_next_round:
                set_of_sg_to_use = set_of_sg_to_use.union(sgrna_for_next_round[sorted_subset])
    else:
        set_of_sg_to_use = set()
        subgroups_sgs = alternative[0]
        candidate = alternative[1]
        for subset in subsets:
            sorted_subset = tuple(sorted(list(subset)))
            if sorted_subset in subgroups_sgs:
                sgs_to_add = subgroups_sgs[sorted_subset]
                for sg in sgs_to_add:

                    if candidate == 1:
                        if sg.seq in sgrna_for_next_round:
                            set_of_sg_to_use.add(sg.seq)
                    else:
                        if sg in sgrna_for_next_round:
                            set_of_sg_to_use.add(sg)


    list_of_sg_to_use = list(set_of_sg_to_use)
    #list_of_sg_to_use.sort(reverse = True, key = lambda x: dict_of_threshold_scores[x])
    list_of_sg_to_use = sort_according_to_score(list_of_sg_to_use, dict_scores, dict_of_threshold_scores) #TODO: isn't it redundant?
    return list_of_sg_to_use
##########################################################################################

def sort_according_to_score(list_of_sg, dict_scores, dict_of_threshold_scores):
    final_list_of_genes = []
    all_scores_sgs_dict = {}
    for sg in list_of_sg:
        score = dict_of_threshold_scores[sg]
        if score in all_scores_sgs_dict:
            all_scores_sgs_dict[score].add(sg)
        else:
            all_scores_sgs_dict[score] = {sg}
    sorted_scores = sorted(list(all_scores_sgs_dict.keys()), reverse = True)
    for score in sorted_scores:
        list_of_sgs_per_score = list(all_scores_sgs_dict[score])
        list_of_sgs_per_score.sort(key = lambda x: dict_scores[x], reverse= True)
        final_list_of_genes.extend(list_of_sgs_per_score)
    return final_list_of_genes
########################################################################################################################
def add_seq_to_gene_set(subgroups_and_sgs, gene_set, candidate):
    if candidate == 1:
        return set([cand.seq for cand in subgroups_and_sgs[gene_set]])
    return set(subgroups_and_sgs[gene_set])


########################################################################################################################
def create_dict_with_dummy_values(set_of_all_possible_subgroups):
    dic = {}
    for subset in set_of_all_possible_subgroups:
        dic[subset] = None
    return dic

########################################################################################################################
# new version of the function. The function filters sgRNAs with a less strict condition
# Param1: dictionary- key = sgRNA seq, value = candidate instance
# Param2: Crispys threshold
# Param3: dictionary of where key = sgRNA sequence, value = (set of genes, dictionary of positions, overall score)
# Param4: dictionary- key = gene, value = list of exons
# Param5: threshold dstance for overlap
# Param6: set of removed sgRNAs
# Param7: sgRNA1 seq
# Param8: sgRNA2 seq
# Param9: List of still remained candidates
# No return: updated sgRNA list (after filtering)
def decide_whether_to_filter(dict_of_sg_candidate, crispys_threshold, dict_of_pos,
                             dict_of_cds, threshold_distance, removed_sgs, sgrna1, sgrna2, proccessed_sg_list):

    candidate1 = dict_of_sg_candidate[sgrna1]
    candidate2 = dict_of_sg_candidate[sgrna2]
    genes_and_scores_sg1 = candidate1.genes_score_dict
    genes_and_scores_sg2 = candidate2.genes_score_dict

    genes_with_low_score = []
    count_worse_score = 0


    for gene in genes_and_scores_sg2:

        if not gene in genes_and_scores_sg1:
            gene_score_of_sg1 = 0
        else:
            gene_score_of_sg1 = genes_and_scores_sg1[gene]
        gene_score_of_sg2 = genes_and_scores_sg2[gene]

        if gene_score_of_sg1 < crispys_threshold and gene_score_of_sg2 < crispys_threshold:
            if gene_score_of_sg2 <= gene_score_of_sg1:
                genes_with_low_score.append(1) #worse
            else:
                genes_with_low_score.append(0)

        elif gene_score_of_sg2 < gene_score_of_sg1:
            count_worse_score += 1

    number_of_high_score_genes = len([gene for gene in genes_and_scores_sg2 if genes_and_scores_sg2[gene] >= crispys_threshold])

    worse_scored = check_if_worse_score(count_worse_score, number_of_high_score_genes, genes_with_low_score, candidate1, candidate2, crispys_threshold)

    overlapping = check_if_overlapping_positions(sgrna1, sgrna2, dict_of_pos, dict_of_cds, threshold_distance, crispys_threshold, genes_and_scores_sg2)
    if worse_scored and overlapping:
        #print("sgRNA is removed ?????????????????????????????  " + str(candidate2.seq) + "\n")
        proccessed_sg_list.remove(candidate2)
        removed_sgs.add(candidate2.seq)


########################################################################################################################
# Param1: number of genes where there is a worse score for sg2 (when sg1 is above the threshold)
# Param2: number of genes with score above threshold for sg2
# Param3: a list for scores beneath the threshold for sg1 ann sg2. 1 if the score for sg2 is worse. Otherwise- 0.
# Param4: sg candidate 1
# Param5: sg candidate 2
# Param6: Crispys threshold (omega)
# Return value- True if sgRNA2 (candidate2) is a candidate for filtering. Otherwise- False.

def check_if_worse_score(count_worse_score, number_of_high_score_genes, genes_with_low_score,
                         candidate1, candidate2, crispys_threshold):
    worse_scored = False
    if count_worse_score > THRESHOLD_FRAC_WORSE_SCORE * number_of_high_score_genes:
        worse_scored = True
    elif count_worse_score == THRESHOLD_FRAC_WORSE_SCORE * number_of_high_score_genes:
        scores_above_thr_sg1 = np.array([candidate1.genes_score_dict[gene] for gene in candidate1.genes_score_dict \
                                         if candidate1.genes_score_dict[gene] >= crispys_threshold])
        scores_above_thr_sg2 = np.array([candidate2.genes_score_dict[gene] for gene in candidate2.genes_score_dict\
                                         if candidate2.genes_score_dict[gene] >= crispys_threshold])
        if np.sum(scores_above_thr_sg1) > np.sum(scores_above_thr_sg2):
            worse_scored = True
            return worse_scored

        elif np.sum(np.array(genes_with_low_score)) > THRESHOLD_FRAC_WORSE_SCORE * len(genes_with_low_score):
            worse_scored = True
    else:
        worse_scored = False
    return worse_scored


########################################################################################################################
def create_copy_set_of_sgs(sets_of_sgRNAs, removed_sgs, copy_set_of_sgRNAs):
    for key in sets_of_sgRNAs:

        sgs = sets_of_sgRNAs[key]
        for sg in sgs:
            if not sg in removed_sgs:
                if not key in copy_set_of_sgRNAs:
                    copy_set_of_sgRNAs[key] = {sg}
                else:
                    copy_set_of_sgRNAs[key].add(sg)

    return



#############################################################################################################
# new version of the filtering function. This function is more strict- will filter much more sgRNAs
# Param1: sgRNA seq #1
# Param2: sgRNA seq #2
# Param3: dictionary of where key = sgRNA sequence, value = (set of genes, dictionary of positions, overall score)
# Param5: dictionary where key = gene, value = exons sequences
# Param5: threshold distance (float)
# Param6: Crispys threshold (omega)
# Param7: dictionary of genes with scores
# Return value: True if all the genes' positions are considered to overlap.
def check_if_overlapping_positions(sgRNA1, sgRNA2, dic, cds_dic, threshold_distance, crispys_threshold, genes_score_sg2):

    positions_sgRNA1 = dic[sgRNA1][1]
    positions_sgRNA2 = dic[sgRNA2][1]
    threshold_distance = int((1-threshold_distance)*20)
    counter = 0
    only_high_scored_genes = [gene for gene in positions_sgRNA2 if genes_score_sg2[gene] >= crispys_threshold]
    number_of_genes_above_thr = len(only_high_scored_genes)    # genes

    for gene in only_high_scored_genes:   #all genes of sgrna2 are contained in sgrna1
        counter_per_gene = 0    # counts overlapped positions
        for pos2 in positions_sgRNA2[gene]:
            #if not gene in positions_sgRNA1:    ##if the score is zero and the format is not the same. example:AT1G13090 (0)
                #return False
            if not gene in positions_sgRNA1:    # don't give penalty
                continue

            for pos1 in positions_sgRNA1[gene]:
                if check_positions_for_overlap(pos2, pos1, threshold_distance, cds_dic[gene]):
                    counter_per_gene += 1
                    break
        # if 90% of the gene's positions overlap, we will considered it as an overlap for a given gene.
        if counter_per_gene > THRESHOLD_NUM_SAME_POSITONS_PER_GENE * len(positions_sgRNA2[gene]):
            counter += 1

    if counter == number_of_genes_above_thr:    # return true if all the genes' positions overlap. Else return False.
        return True
    return False

##########################################################################################################


def remove_sgs_from_list_of_candidates(sgs_per_internal_node_dic,  first_candidates):
    new_lst = []
    for gene_set in sgs_per_internal_node_dic:

        sgrnas = sgs_per_internal_node_dic[gene_set]
        for sgrna in sgrnas:
            for candidate in first_candidates:
                if candidate.seq == sgrna:
                    new_lst.append(candidate)
                    break
    return new_lst

######################################################################################################

#A function which creates a set on- target genes which the sgRNA cuts according to the output file.
#Input: genes and all the related data (including target genes, score and position).
#Output: set of on- target genes.


def extract_genes(genes_and_pos):

    #p = re.compile("([\S]+)[\s]+\(")

    lst_of_genes = []
    lst_of_positions_per_gene = []
    p = re.compile("([\S]+)[\s]+\((.*?)\)")

    pattern = re.compile("pos:[\s]+([\d]+R{0,1})")

    lst_of_genes_positions = p.findall(genes_and_pos)
    for gene_pos in lst_of_genes_positions:
        lst_of_genes.append(gene_pos[0])
        lst_of_positions_per_gene.append(pattern.findall(gene_pos[1]))

    genes_set = set(lst_of_genes)
    dict_of_genes_and_pos = {}

    for i in range(len(lst_of_genes)):
        for position in lst_of_positions_per_gene[i]:

            if position[-1] == 'R':
                if not lst_of_genes[i] in dict_of_genes_and_pos:
                    dict_of_genes_and_pos[lst_of_genes[i]] = [(int(position[:-1]), False)]
                else:
                    dict_of_genes_and_pos[lst_of_genes[i]].append((int(position[:-1]), False))

            else:

                if not lst_of_genes[i] in dict_of_genes_and_pos:
                    dict_of_genes_and_pos[lst_of_genes[i]] = [(int(position), True)]
                else:
                    dict_of_genes_and_pos[lst_of_genes[i]].append((int(position), True))

    return genes_set, dict_of_genes_and_pos

#######################################################################################################################

# Param1: postions of sgRNA1
# Param2: positions if sgRNA2
# Param3: threshold overlapping distance
# Param4: cds exons
# Return Value: 1 if the positions are overlapped, and 0 if not.

def check_positions_for_overlap(pos_sgRNA1, pos_sgRNA2, threshold_distance, cds_exons):
    #distance = 0

    #introns_lengths = get_introns_lengths(exons)
    strand_sgRNA1 = pos_sgRNA1[1]
    strand_sgRNA2 = pos_sgRNA2[1]
    start_sgRNA1 = pos_sgRNA1[0]
    start_sgRNA2 = pos_sgRNA2[0]

    cds = ""
    for exon in cds_exons:
        cds += exon
    length_cds = len(cds)

    if strand_sgRNA1 == strand_sgRNA2:
        if strand_sgRNA1:
            new_pos_sgRNA1 = start_sgRNA1
            new_pos_sgRNA2 = start_sgRNA2
        else:
            new_pos_sgRNA1 = len(cds)-(start_sgRNA1+20) #reverse indices
            new_pos_sgRNA2 = len(cds)-(start_sgRNA2+20)

    else:

        if strand_sgRNA1:

            new_pos_sgRNA1 = start_sgRNA1
            new_pos_sgRNA2 = length_cds-(start_sgRNA2+20)
        else:
            new_pos_sgRNA1 = length_cds-(start_sgRNA1+20)
            new_pos_sgRNA2 = start_sgRNA2

    different_exons = check_if_different(new_pos_sgRNA1, new_pos_sgRNA2, cds_exons)
    # first check if the positions are on the same exon. The cds consits of concatenated exons, so it could be
    # that the positions are closed relative to the cds sequence, but pratically they could be distant (on different
    # exons).
    if different_exons: # no overlap
        return 0

    # There is an overlap
    distance = abs(new_pos_sgRNA1-new_pos_sgRNA2)   #The bigger is the distance, the smaller region is overlapped
    # suppose the parameter of the overlap threshold in filter_sgRNA_main.py was 0.1, it means that
    # the threshold distance is 18 nuc (1-threshold) *20).
    # If the distance between the positions is less than 18, it means that there is an overlap of more than 0.1 (2 bp)
    # betweeen the two examined positions.
    if distance <= threshold_distance:  #distance is less or equal to the threshold-> large overlapped region
        return 1    #penalty
    return 0
###################################################################################################
# Param1: position in cds for sgRNA1
# Param2: position in cds for sgRNA2
# Param3: list of exons for a given gene
# Return value: True if positions are on different exons. Else: False.
def check_if_different(new_pos_sgRNA1, new_pos_sgRNA2, cds_exons):
    exon_sgrna_1 = 0
    exon_sgrna_2 = 0
    lst_of_indices = []
    for i in range(len(cds_exons)):
        if i == 0:
            lst_of_indices.append((0, len(cds_exons[0])-1))
        else:
            lst_of_indices.append((lst_of_indices[i-1][1]+1, lst_of_indices[i-1][1]+1+len(cds_exons[i])-1))
        if new_pos_sgRNA1 >= lst_of_indices[i][0] and new_pos_sgRNA1 <= lst_of_indices[i][1]:
            exon_sgrna_1 = i
        if new_pos_sgRNA2 >= lst_of_indices[i][0] and new_pos_sgRNA2 <= lst_of_indices[i][1]:
            exon_sgrna_2 = i

    if exon_sgrna_1 != exon_sgrna_2:
        return True
    return False

###########################################################################################################

# Returns two items:
# 1. a dictionary with a tuple of genes as a key and a sorted list of sets of genes which are subsets of the key tuple
# 2. a sorted list of sets (sorted according to length)
def create_dict_of_subsets(sets_of_sgRNAs, bottom_up = False, top_down = False):
    dic = {}
    lst_of_sets = list(sets_of_sgRNAs.keys())
    #sorted_list = sorted(lst_of_sets, reverse = not bottom_up)
    sorted_list = sorted(lst_of_sets, key = lambda x: len(x), reverse = not bottom_up)  # reverse will be True if the parameter is top down.
    for i in range(len(sorted_list)):
        for j in range(len(sorted_list)):
            if (set(sorted_list[j]) <= set(sorted_list[i]) and (not top_down))\
                    or ((set(sorted_list[j]) >= set(sorted_list[i])) and top_down) :  # always need a list of subsets
                if not sorted_list[i] in dic:
                    dic[sorted_list[i]] = [set(sorted_list[j])]
                else:
                    dic[sorted_list[i]].append(set(sorted_list[j]))
        #dic[sorted_list[i]].sort(reverse = not bottom_up)
    #sorted_list.sort(reverse= not bottom_up, key = lambda x: len(x))
    return dic, sorted_list
#########################################################################################################3##############
def create_dict_from_sg_list(list_of_sg):
    dic = {}
    for i in range(len(list_of_sg)):
        sg = list_of_sg[i].seq
        dic[sg] = i
    return dic
#######################################################################################################################


def leave_only_relevant_sgRNA(res):
    if len(res) < 1:
        return
    #candidates_to_del = set()
    for i in range(len(res) -1, -1, -1):
        if res[i].cut_expectation < 1:

            #candidates_to_del.add(res[i])
            del res[i]

        elif i < len(res) -1:
            for j in range(i+1, len(res)):
                if j >= len(res):
                    continue
                if res[i].seq == res[j].seq:
                    if res[i].cut_expectation <= res[j].cut_expectation:

                        #candidates_to_del.add(res[i])
                        del res[i]
                    else:
                        #candidates_to_del.add(res[j])
                        del res[j]

    return #candidates_to_del
#######################################################################################################################

def delete_singletons(proccessed_sg_list, threshold):
    copied_list = proccessed_sg_list.copy()
    for candidate in copied_list:
        scores_of_genes = list(candidate.genes_score_dict.values())
        remain = 0
        for score in scores_of_genes:
            if score >= threshold:
                remain += 1
        if remain < 2:
            #print("**************************  REMOVE SINGLETON *********************** "+str(candidate.seq)+"\n")
            proccessed_sg_list.remove(candidate)
    return


