from create_dictionaries_and_pickle import *
from sgRNA_classification import *
import argparse
import time
import shutil
from global_variables import *

# Param1: full path to the family directory (or to the directory of the second round run)
# Param2: family name
# Param3: annotation path for the family genes
# Param4: a path for the dictionary of chromosomes and their corresponding lengths.
# Param5: full path to the family directory.
# Param6: threshold for the off-target score.
# Param7: threshold distance- threshold for the proportion of the overlap between two sgRNAs.
# Param8: The required number of sgRNAs. If is zero, there is no limitation.
# param9: 1 if bottom up. 0 otherwise
# Param10: 1 if there is a need for a second run. 0 otherwise.
# Param11: Cripys threshold (omega).
# Param12: If the threshold for off-targets is relative, this parameter is set to 1. Otherwise, if it is absolute,
# it is set 0.

def filter_sgRNAs_family(full_dir_path, family, annotation_file_path, dic_genome_path, direct_family_dir, threshold,
                         threshold_distance, number_of_sgs, bottom_up, method, crispys_threshold, relative):
    if bottom_up == 1:
        bottom_up = True
    else:
        bottom_up = False
    #set_of_used_sg = set()
    low_threshold = LOW_THRESHOLD_OMEGA  # can be changed later if we want to make just one run.
    if method == 0: # if I want to run the filtering procedure only once
        low_threshold = crispys_threshold

    start_time = time.time()    #measure running time

    lst_of_all_sg_pickle = os.path.join(full_dir_path, 'res_in_lst.p')
    lst_sgRNAs = load_pickle_files(lst_of_all_sg_pickle) #A list of lists of sgRNAs (divided into sets)
    if lst_sgRNAs == []:
        return
    set_of_subgroups = create_set_of_internal_nodes(lst_sgRNAs) # finds the true internal nodes according to the output
    # proccessed_sg_list:
    # list where each element is: instance of class Candidate:
    # candidate.seq = sgRNA sequence
    # candidate.cut_expectation = the overall score of the sgRNA (sum of all genes' score)
    # candidate.genes_score_dict = dictionary, where key = gene, value = score per gene
    # the last printed doctionary: key = gene, value = a list of lists, where each element is
    # [sgRNA seq, overall score, dictionary of genes and score, dictionary with genes and aligned region with substitutions]
    # specifically, 'dictionary of genes and score' is: key = gene, value = score per gene.
    # 'dictionary with genes and aligned region with substitutions': key = gene, value = [aligned region, dictionary
    # (key = index of mm, value = [char in sgRNA, char in the aligned region])]
    proccessed_sg_list = concatenate_list(lst_sgRNAs)    # concatenates all the subsets
    initial_number_of_sgs_before_del = len(set([candidate.seq for candidate in proccessed_sg_list]))
    delete_singletons(proccessed_sg_list, crispys_threshold)  # eliminating sgrnas which bind to a single gene--> need sgrnas for multiple genes
    leave_only_relevant_sgRNA(proccessed_sg_list)    #first filtering of sgRNAs (removing repetitions)
    initial_number_of_sgs = len(proccessed_sg_list) #afer the deletions of single targeted sgRNAs
    dict_genes_score_for_sgrnas = create_dict_scores_per_sg(proccessed_sg_list) #scores of all the genes (not only above threshold)
    dic_genome = load_pickle_files(dic_genome_path) # a dictionary of chromosomes and sequence lengths

    #family_new_dir_path = os.path.join(direct_family_dir, family)
    cds_file = os.path.join(full_dir_path, family + ".txt") # a path for the gene file
    dict_of_cds = extract_cds(cds_file) # key is gene name and value is a list of exons

    file_path = os.path.join(full_dir_path, "CRISPys_output.csv")
    file = open(file_path, 'r')
    content = file.read()
    file.close()

    dict_sg_candidate = create_dict_sg_candidate(proccessed_sg_list) # dict where key = sgRNA, value = candidate instance
    modified_content, entries = extract_sgrna_info(content, dict_sg_candidate)  # here is the updated function which parses the file
    list_of_sgRNAs_remained = []
    dictionary_of_remained_lines = {}
    count_non_filtered_off_target = 0
    #indicator_CRISTA = os.path.exists(os.path.join(full_dir_path, 'CRISTA'))
    sets_of_sgRNAs, dict_of_sg_genes_pos, dict_sgrna_score_threshold = {}, {}, {}
    family_new_dir_path_sgRNA = os.path.join(direct_family_dir, DIR_NAME_BWA, "off_target_results")
    if not os.path.exists(family_new_dir_path_sgRNA):
        raise Exception("filter_sgRNAs_family() in filter_sgRNAs.py: BWA was not run for this family\n")

    ## starting the filtering procesure of off-target and position overlap
    for sgRNA in modified_content:
        score = modified_content[sgRNA][0]
        genes_and_positions = extract_correct_orientations(modified_content[sgRNA][1])
        genes = set(list(genes_and_positions.keys()))
        #sgRNA_data = modified_content[sgRNA][2]
        sgRNA_data = modified_content[sgRNA][3]
        #family_new_dir_path_sgRNA = os.path.join(direct_family_dir, DIR_NAME_BWA, "sgRNA_" + sgRNA)
        file_path_BWA_res = find_crista_output_filename(family_new_dir_path_sgRNA, sgRNA)
        compatible_sgRNA = False
        if not file_path_BWA_res:
            compatible_sgRNA = True
        if not compatible_sgRNA:
            dict_of_BWA_results = parse_crista_results(file_path_BWA_res) # get a dictionary of BWA results
            compatible_sgRNA = check_sgRNAs_compatability_alternative(threshold,
                                                                      relative, dict_genes_score_for_sgrnas[sgRNA],
                                                                      dict_of_BWA_results, annotation_file_path,
                                                                       dic_genome)  # , path_for_info)
        if compatible_sgRNA:
            # If the sgRNA is compatible several things occur:
            # 1. Update the counter of non-filtered sgRNAs
            # 2. Update the list of remained sgRNAs
            # 3. Update the lines to be written to the csv result file.
            # 4. Clustering the new sgRNA to the relevant subgroup of genes (the ones with a score above threshold).
            # 5. Update the dictionary where key= sgRNA seq, value = the total score (considering only the genes>= threshold
            candidate_new_sg = dict_sg_candidate[sgRNA] #find_candidate(sgRNA, proccessed_sg_list)
            count_non_filtered_off_target += 1
            new_element_to_remained_list = [sgRNA, score, genes, genes_and_positions]
            list_of_sgRNAs_remained.append(new_element_to_remained_list)
            dictionary_of_remained_lines[sgRNA] = sgRNA_data # it is just stored to allow a convienient wrting of the csv results file.
            # update sets_of_sgRNAs and dict_of_sg_genes_pos dictionaries
            cluster_new_sgRNA(sgRNA, new_element_to_remained_list, sets_of_sgRNAs,
                              dict_of_sg_genes_pos, candidate_new_sg, low_threshold)
            update_dict_score_threshold(candidate_new_sg, crispys_threshold, dict_sgrna_score_threshold)

        else:
            print("Not compatibale sgRNA "+ sgRNA+"\n") # a print for log
            proccessed_sg_list.remove(dict_sg_candidate[sgRNA])    #remove incmpatible sg
    print("after off-target search:")
    print(proccessed_sg_list)

    # filter sgRNAs with overlapping positions.
    subgroups_and_sgs = filter_same_positioned_sg(sets_of_sgRNAs, dict_of_cds, proccessed_sg_list,
                                                  dict_of_sg_genes_pos, threshold_distance, crispys_threshold) #updated set of clustered sgRNAs per subgroup
    not_filtered_overlaping = len(proccessed_sg_list)
    print("after filtering same positioned:")
    print(proccessed_sg_list)
    #print("**************************** After deleting overlaps length is " + str(len(proccessed_sg_list)) + "\n")
    final_number_of_sgs = None
    dic_scores = create_dict_score_sgs(proccessed_sg_list)  # a dictionary- key = sgRNA seq, value = the total score
    if number_of_sgs != 0:  # if we have a particular limitation of sgRNA per internal node
        sgs_per_internal_node_dic = create_dict_of_sgs_per_subgroup_alternative(subgroups_and_sgs, 0, number_of_sgs,
                                                                                dic_scores, set_of_subgroups,
                                                                                bottom_up,
                                                                                dict_sgrna_score_threshold)
        less_than_needed_sg = [gene_set for gene_set in sgs_per_internal_node_dic if len(
            sgs_per_internal_node_dic[gene_set]) < number_of_sgs]
        num_internal_nodes_in_dic = len(extract_only_internal_nodes(sgs_per_internal_node_dic, set_of_subgroups))
        number_of_internal_nodes = len(dict_of_cds) - 1
        if (not "second_run" in full_dir_path) and method:
            if (num_internal_nodes_in_dic < number_of_internal_nodes or len(less_than_needed_sg) > 0) \
                    and (crispys_threshold != low_threshold):
                new_output_directory = os.path.join(full_dir_path, 'second_run')
                copy_second_round_crispys(new_output_directory, full_dir_path, family)

            else:  # No need to run second round
                ## seems like this function is for another sanity check...
                final_sgRNAs_list = remove_sgs_from_list_of_candidates(sgs_per_internal_node_dic, proccessed_sg_list)
                final_number_of_sgs = len(final_sgRNAs_list)
        path_for_internal_node_dict = os.path.join(full_dir_path, 'internal_nodes_per_sg.p')
        pickle.dump(sgs_per_internal_node_dic, open(path_for_internal_node_dict, 'wb'))
    else:
        final_sgRNAs_list = proccessed_sg_list
        final_number_of_sgs = len(final_sgRNAs_list)


    # remove the filtered sgRNAs from the list of remained sgRNAs and dictionary of the remained lines.
    remove_unneeded_sg(list_of_sgRNAs_remained, proccessed_sg_list, dictionary_of_remained_lines)
    path_for_filtered_statistics = os.path.join(full_dir_path, 'number_of_sgs.p')
    dict_of_numbers_of_sgs = {0 : initial_number_of_sgs_before_del, 1: initial_number_of_sgs,
                              2: count_non_filtered_off_target, 3: not_filtered_overlaping}

    if final_number_of_sgs != None:

        dict_of_numbers_of_sgs[4] = final_number_of_sgs
        path_for_final_list = os.path.join(full_dir_path, 'res_in_filtered_lst_final.p')
        path_for_final_csv_file = os.path.join(full_dir_path, 'output_filtered_final.csv')
        pickle.dump(final_sgRNAs_list, open(path_for_final_list, 'wb'))

    pickle.dump(dict_of_numbers_of_sgs, open(path_for_filtered_statistics, 'wb'))
    # write results
    write_output_file(os.path.join(full_dir_path, "output_filtered.csv"), list_of_sgRNAs_remained,
                      entries, dictionary_of_remained_lines)
    current_time = time.time()
    run_time = current_time-start_time
    print("runtime is : "+str(run_time)+" sec\n")
    print("The job has finished correctly!")
    path_for_pickle_remained_list = os.path.join(full_dir_path, "res_in_filtered_lst.p")
    pickle.dump(proccessed_sg_list, open(path_for_pickle_remained_list, 'wb'))
    if final_number_of_sgs != None:
        remove_unneeded_sg(list_of_sgRNAs_remained, final_sgRNAs_list, dictionary_of_remained_lines)
        write_output_file(path_for_final_csv_file, list_of_sgRNAs_remained, entries,
                          dictionary_of_remained_lines)

    return

######################################################################################################################
def create_dict_of_sgrna_genes_above_threshold(proccessed_sg_list, crispys_threshold):
    dic_sg_genes_threshold = {}
    for candidate in proccessed_sg_list:
        dict_gene_score = candidate.genes_score_dict
        dic_sg_genes_threshold[candidate.seq] = set()
        for gene in dict_gene_score:
            if dict_gene_score[gene] >= crispys_threshold:
                dic_sg_genes_threshold[candidate.seq].add(gene)
    return dic_sg_genes_threshold


#######################################################################################################################
# A function which get a list of elements of type Candidate as input
# Returns: a dictionary: key = sgRNA sequence, value = list of scores.
def create_dict_scores_per_sg(proccessed_sg_list):
    dic = {}
    for candidate in proccessed_sg_list:
        dic[candidate.seq] = list(candidate.genes_score_dict.values())
    return dic
########################################################################################################################
# update a dictionary dict_sgrna_score_threshold in the following manner: key = sgRNA seq, value= overall score of
# sgRNA, composed only of genes above threshold
def update_dict_score_threshold(candidate, threshold, dict_sgrna_score_threshold):
    score_of_sg = 0
    dict_gene_score = candidate.genes_score_dict
    for gene in dict_gene_score:
        if dict_gene_score[gene] >= threshold:
            score_of_sg += dict_gene_score[gene]
    dict_sgrna_score_threshold[candidate.seq] = score_of_sg
    return
#######################################################################################################################
def delete_sg_from_everywhere(candidate, proccessed_sg_list, set_of_sgRNAs_to_fill,
                              dictionary_of_remained_lines, sets_of_sgRNAs, dict_of_sg_genes_pos,
                              dict_sgrna_score_threshold):
    proccessed_sg_list.remove(candidate)
    set_of_sgRNAs_to_fill -= {candidate}
    dictionary_of_remained_lines.pop(candidate.seq)
    gene_set = tuple(sorted(list(candidate.genes_score_dict.keys())))
    sgRNAs_of_set = sets_of_sgRNAs[gene_set]
    sgRNAs_of_set -= {candidate.seq}
    dict_of_sg_genes_pos.pop(candidate.seq)
    dict_sgrna_score_threshold.pop(dict_sgrna_score_threshold)
    return

#######################################################################################################################

def remove_unneeded_sg(list_of_sgRNAs_remained, proccessed_sg_list, dictionary_of_remained_lines):


    copied_list_of_remained = list_of_sgRNAs_remained.copy()
    for item in copied_list_of_remained:
        sg = item[0]
        bool = False
        for candidate in proccessed_sg_list:
            if candidate.seq == sg:
                bool = True
                break
        if not bool:
            list_of_sgRNAs_remained.remove(item)
            sg_del = dictionary_of_remained_lines.pop(sg)
            #print("sgRNA was removed : "+sg_del+"!\n")

########################################################################################################################
# This function updated the following data structures:
# 1. sets_of_sgRNAs: a dictionary- key = sorted subgroup of genes, value = a set of the corresponding sgRNAs.
# 2. dict_of_sg_genes_pos: a dictionary- key = sgRNA sequence, value = (set of genes, dictionary of positions, overall score)
#   where dictionary of positions is: key = gene, value = list of positions where each element is a tuple (position, boolean),
#   The boolean is True if the strand is + and Flase if the strand is -.
# Input1: sgRNA sequence
# Input2: a list [sgRNA seq, overall score, set of genes, dictionary of genes (keys) and positions (value) represented
#   by tuples (position, boolean) where boolean is True when the strabd is + and False when the strand is -.
# In sets_of_genes the key is a sorted tuple of genes that are above (or equal) the threshold score.

def cluster_new_sgRNA(sgRNA, new_element_to_remained_list, sets_of_sgRNAs, dict_of_sg_genes_pos,
                      candidate_new_sg, low_threshold):

    score = new_element_to_remained_list[1]
    genes = new_element_to_remained_list[2]
    positions = new_element_to_remained_list[3]

    dict_of_sg_genes_pos[sgRNA] = (genes, positions, float(score))
    genes_above_threshold = [gene for gene in candidate_new_sg.genes_score_dict \
                             if candidate_new_sg.genes_score_dict[gene] >= low_threshold]
    sorted_genes = tuple(sorted(genes_above_threshold))

    if not sorted_genes in sets_of_sgRNAs:
        sets_of_sgRNAs[sorted_genes] = {candidate_new_sg.seq}
    else:
        sets_of_sgRNAs[sorted_genes].add(candidate_new_sg.seq)
    return

#######################################################################################################################

# get a list of sgRNAs (altogether)
def concatenate_list(lst_sgRNAs):
    concat_list = []    #concatenated list of all the subsets
    for lst in lst_sgRNAs:
        concat_list.extend(lst.candidate_lst)
    return concat_list




#######################################################################################################################

def extract_correct_orientations(dict_genes_and_positions):
    dict_genes_positions_orientation = {}
    labels = {'+':True, '-': False}
    for gene in dict_genes_and_positions:
        dict_genes_positions_orientation[gene] = []
        positions = dict_genes_and_positions[gene]
        for position in positions:
            dict_genes_positions_orientation[gene].append((int(position[:-1]), labels[position[-1]]))
    return dict_genes_positions_orientation

########################################################################################################################33

def check_if_remained_in_list(sgRNA, genes, proccessed_sg_list):
    possible_candidates = []
    remained = False

    for candidate in proccessed_sg_list:
        if candidate.seq == sgRNA:
            list_genes = candidate.genes_score_dict.keys()
            non_zero_scored_genes = set([gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] > 0])
            possible_candidates.append(non_zero_scored_genes)
            set_of_cand_genes = set(list(list_genes))
            #possible_candidates.append(set_of_cand_genes)
            if set_of_cand_genes == genes:

                remained = True #to remain
                break
    if not remained:
        #if sgRNA == "AGGAGAAAGAGAGAAGGAGA":
            #print("--------- set of candidates is   ------------\n")
            #print("genes are : "+ str(genes)+"\n")
            #print("set of candidates are : "+"\n")
        for possible_candidate in possible_candidates:
            #if sgRNA == "AGGAGAAAGAGAGAAGGAGA":
                #print(str(possible_candidate)+"\n")

            if possible_candidate == genes:
                remained = True
                break


    return remained



########################################################################################################
# Input:
# Param1: output file path (csv file)
# Param2: list of sgRNAs to write
# Param3: Titles for the table in the output file.
# Output:
# Return: No return value. The function updates the output file passed by the first argument.

def write_output_file(output_file_path, list_of_sgRNAs, entries, dict_remained_lines):
    list_of_sgRNAs.sort(key = lambda x: x[1], reverse=True)

    file = open(output_file_path, 'w')
    #file.write(entries+"\n")
    file.write(entries+"\n")
    for data in list_of_sgRNAs:
        sg = data[0]
        if sg in dict_remained_lines:
            file.write(dict_remained_lines[sg])
            file.write("\n")

    file.close()




################################################################################


def update_sg_subsets_dict(sg_subgroup_dic, sgrnas, gene_set):
    for sg in sgrnas:

        sg_subgroup_dic[sg] = gene_set


########################################################################################################################
# Param1: the dst directory where the new Crispys run should be executed.
# Param2: the directory of the first run- the exons fasta file is copied from there.
# param3: family name
# No return value.
# In case the global variable of the second run was specified, the results files will be coppied from there (but not the
# BWA directory).
def copy_second_round_crispys(new_output_directory, main_family_dir, family):

    if os.path.exists(new_output_directory):
        return  # means I have already created the folder and copied all the relevant files

    os.makedirs(new_output_directory)
    second_round_dst_folder = os.path.join(new_output_directory, family)
    os.makedirs(second_round_dst_folder)
    shutil.copyfile(os.path.join(main_family_dir, family + ".txt"),
                    os.path.join(second_round_dst_folder, family + ".txt"))
    if os.path.exists(os.path.join(main_family_dir, "protdist")):
        shutil.copyfile(os.path.join(main_family_dir, "protdist"), os.path.join(second_round_dst_folder, "protdist"))

    return

########################################################################################################################

# Among all the found genes sets of candidate find the ones which really represent internal nodes

def extract_only_internal_nodes(sgs_per_internal_node_dic, set_of_subgroups):
    subsets_genes = list(sgs_per_internal_node_dic.keys())
    internal_nodes = [internal_node for internal_node in subsets_genes if internal_node in set_of_subgroups]
    return internal_nodes


########################################################################################################################
# A function that combines information from both the dictionary of candidates and the content from csv Crispys output file
# Note: The csv contains also the positions of the sgRNA alignment with the target sequences.
# Input1: content from Crispys csv file (string)
# Input2: dictionary, where key = sgRNA seq, value = Candidate instance
# Return val1: a dictionary, where key = sgRNA seq, value = list of the following elements
#   1. Element 1: the overall score
#   2. Element 2: dictionary, where key = gene, value = list of positions (relative to the cds length). The positions
#       can be on the forward strand (+) or on the reverse strand (-).
#   3. Element 3: dictionary, where key = gene, value = score
#   4. Element 4: the entry as it appears in the csv file (string)
# Return val2: The header of the Crispys output csv file.

def extract_sgrna_info(content, dict_sg_candidate):
    number_of_parameters = 8
    header = None
    dict_of_sg_and_info = {}
    # splitting to tables- in Crispys each table id represents a different subset of genes
    tables = re.split('table id', content)
    for i in range(len(tables)):
        if i == 0:
            continue
        rows = re.split('\n', tables[i])
        for j in range(len(rows)):
            if rows[j] == '':
                continue
            elif j == 1 or j == 0:
                continue
            splitted_row = re.split(",", rows[j])
            if len(splitted_row) != number_of_parameters:
                continue
            if not header:
                header = rows[1]
            if splitted_row[1] != '':
                sgrna = splitted_row[1]
                score = splitted_row[2]
            if splitted_row[3] != '':
                gene = splitted_row[3]
            pos = splitted_row[7]
            #print("*"+splitted_row[4]+"*\n")
            gene_score = splitted_row[4]
            if gene_score != "":
                gene_score = float(splitted_row[4])
            if not (sgrna, i) in dict_of_sg_and_info:
                dict_of_sg_and_info[(sgrna, i)] = [score, {gene: [pos]}, {gene: gene_score}, rows[j]+"\n"]
            else:
                genes_pos_dict = dict_of_sg_and_info[(sgrna, i)][1]
                genes_score_dict = dict_of_sg_and_info[(sgrna, i)][2]
                dict_of_sg_and_info[(sgrna, i)][3] += rows[j]+"\n"
                if gene_score != "":
                    genes_score_dict[gene] = gene_score
                if gene in genes_pos_dict:
                    genes_pos_dict[gene].append(pos)

                else:
                    genes_pos_dict[gene] = [pos]
    dict_of_sg_and_info_final = edit_dict_of_sgrna_data(dict_of_sg_and_info, dict_sg_candidate)
    return dict_of_sg_and_info_final, header

#######################################################################################################################
# Input1: dictionary, where key = (sgRNA seq, table id), value= [score, {gene: gene score}, the entry]
# Input2: dictionary, where key = sgRNA seq, value = Candidate instance
# Return value: a dictionary, where key = sgRNA seq, value = list of the following elements
# #   1. Element 1: the overall score
# #   2. Element 2: dictionary, where key = gene, value = list of positions (relative to the cds length). The positions
# #       can be on the forward strand (+) or on the reverse strand (-).
# #   3. Element 3: dictionary, where key = gene, value = score
# #   4. Element 4: the entry as it appears in the csv file (string)

def edit_dict_of_sgrna_data(dict_of_sg_and_info, dict_sg_candidate):
    final_dict = {}
    for sgrna_table in dict_of_sg_and_info:
        sgrna = sgrna_table[0]
        if not sgrna in dict_sg_candidate:
            continue
        candidate = dict_sg_candidate[sgrna]
        genes_of_candidate = set(list(candidate.genes_score_dict.keys()))
        non_zero_candidates_genes = set([gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] > 0])
        sgrna_data = dict_of_sg_and_info[sgrna_table]
        genes_score_dict = sgrna_data[2]
        suggested_genes = set(list(sgrna_data[1].keys()))

        if not sgrna in final_dict:
            if genes_of_candidate == suggested_genes:
                final_dict[sgrna] =  dict_of_sg_and_info[sgrna_table]
            elif non_zero_candidates_genes == suggested_genes:
                final_dict[sgrna] = dict_of_sg_and_info[sgrna_table]
        # else:
        #
        #     #genes_of_candidate = set(list(candidate.genes_score_dict.keys()))
        #     if  suggested_genes == genes_of_candidate or non_zero_candidates_genes == suggested_genes:
        #         final_dict[sgrna] = sgrna_data

    return final_dict

########################################################################################################################
def is_equal(candidate_genes_score_dict, genes_score_dict_csv, non_zero_candidates_genes):
    new_dict = {}
    for gene in non_zero_candidates_genes:
        new_dict[gene] = genes_score_dict_csv[gene]
    if candidate_genes_score_dict == genes_score_dict_csv:
        return True
    return False

########################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filters off")
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--directory", "-d", help = "directory")
    parser.add_argument("--annotation", "-a", help = "annotation file path")
    parser.add_argument("--genome", "-g", help = "genome dictionary file path")
    parser.add_argument("--target", "-t", help = "target directory of families")
    parser.add_argument("--threshold", "-w", help = "threshold for determing if a specific sgRNA is off-target")
    parser.add_argument("--threshold_distance", "-y", help = "threshold for overlapping")
    parser.add_argument("--number_of_sgs", "-n", help = "number of sgRNAs per internal node")
    parser.add_argument("--bottom_up", "-b", help = "specify if the filtering should be bottom up")
    parser.add_argument("--method", "-m", help = "specify if to use the default or the alternative")
    parser.add_argument("--crispys_threshold", "-c", help = "crispys_threshold")
    parser.add_argument("--relative", "-r", help = "relative")

    args = parser.parse_args()
    output = args.output
    family = args.directory
    annotation_file = args.annotation
    genome_dic_path = args.genome
    target_dir = args.target
    threshold = args.threshold
    threshold_dist = args.threshold_distance
    number_of_sgs = args.number_of_sgs
    bottom_up = args.bottom_up
    method = args.method
    crispys_threshold = args.crispys_threshold
    relative = args.relative    # 0 if not relative, 1 if relative

    filter_sgRNAs_family(output, family, annotation_file, genome_dic_path, target_dir, threshold, threshold_dist,
                         int(number_of_sgs), int(bottom_up), int(method), float(crispys_threshold), int(relative))



