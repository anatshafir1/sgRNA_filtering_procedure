#from filter_sgRNAs import *
from sgRNA_classification import *
import argparse

# Param1: the directory of the family (first run)
# Param2: the directory of the second run
# Param3: the name of the family
# Param4: the max number of sgRNAs per internal node
# Param5: threshold distance (for overlapping). This is a float
# Param6: if bottom-up 1, else: 0
# Param7: Crispys threshold (the one used for the second round)
def join_output(family_dir, path_dir_second_run, family, number_of_sgs, threshold_distance, bottom_up, crispys_threshold):
    if bottom_up == 1:
        bottom_up = True
    else:
        bottom_up = False
    cds_file = os.path.join(family_dir, family + ".txt")
    dict_of_cds = extract_cds(cds_file)
    path_for_dic_number_sgs = os.path.join(family_dir, 'number_of_sgs.p')
    with open(path_for_dic_number_sgs, 'rb') as handle:
        dict_number_sgs_for_all = pickle.load(handle)
    handle.close()
    with open(os.path.join(family_dir, 'res_in_filtered_lst.p'), 'rb') as handle:
        first_candidates = pickle.load(handle)
    handle.close()
    with open(os.path.join(path_dir_second_run, family, 'res_in_filtered_lst.p'), 'rb') as handle:
        second_candidates = pickle.load(handle)
    handle.close()
    with open(os.path.join(family_dir, 'res_in_lst.p'), 'rb') as handle:
        lst_of_all_internal_nodes_sgs = pickle.load(handle)
    handle.close()

    #if len(lst_of_all_internal_nodes_sgs) != len(dict_of_cds) - 1:

    # set of internal nodes represented as tuples of genes in a sorted order.
    crispys_internal_nodes = create_set_of_internal_nodes(lst_of_all_internal_nodes_sgs)


    output_first_path = os.path.join(family_dir, 'output_filtered.csv')
    output_second_path = os.path.join(path_dir_second_run, family, 'output_filtered.csv')

    file_first = open(output_first_path, 'r')
    file_second = open(output_second_path, 'r')

    content_first = file_first.read()
    content_second = file_second.read()

    file_first.close()
    file_second.close()

    splitted_first_content = content_first.split('\n')[1:]  # skipping the first line (title)
    splitted_second_content = content_second.split('\n')[1:]     # skipping the first line (title)

    dict_of_lines_first_run = {} #For the first run: a dictionary:
                                    # key = sgRNA seq, value= the corresponding line in the csv file
    dict_of_lines_second_run = {} # For the second run: a dictionary- key = sgRNA, value = the relevant line
                                    # in the csv file
    # updating the two above dictionaries
    fill_dict_lines(splitted_second_content, dict_of_lines_second_run)#, None)
    fill_dict_lines(splitted_first_content, dict_of_lines_first_run)#None)

    # subgroups_and_sgs_first: a dictionary of subgroup as a key (exactly the subset of genes to which the
    # # # # sgRNA bounds) and a list of corresponding candidates as values
    # sgRNAs_first: set of sgRNAs (based on the input list)
    subgroups_and_sgs_first, sgRNAs_first = create_subgoups(first_candidates)
    subgroups_and_sgs_second, sgRNAs_second = create_subgoups(second_candidates)


    # combined list of gRNAs (and all the relevant info) from both the runs
    lst_of_remained_sgs = create_lst_of_remained(dict_of_lines_first_run, dict_of_lines_second_run)
    dict_of_sg_genes_pos = create_pos_dict(lst_of_remained_sgs) # Return value: key = sgRNA seq,
                                                        # value = (genes, {gene: [(position, True/False)}, total score)

    merged_set_of_subgroups = merge_subgroups(subgroups_and_sgs_first, subgroups_and_sgs_second, sgRNAs_first)
    merged_lst_of_candidates = merge_list_of_candidates(first_candidates, second_candidates, sgRNAs_first)



    filter_same_positions_for_joined(merged_set_of_subgroups, dict_of_cds, merged_lst_of_candidates, dict_of_sg_genes_pos,
                                     threshold_distance, crispys_threshold) # make sure sgRNAs of the first group do not overlap with the previous one

    # update all the lists after filtering


    first_candidates = update_candidates_list(sgRNAs_first, merged_lst_of_candidates, None) # is it needed???
    subgroups_and_sgs_first, sgRNAs_first = create_subgoups(first_candidates)
    second_candidates = update_candidates_list(sgRNAs_second, merged_lst_of_candidates, sgRNAs_first)
    subgroups_and_sgs_second, sgRNAs_second = create_subgoups(second_candidates)

    # dictionaries with key = sgRNA seq, value= total score
    dict_of_scores_first = create_dict_score_sgs(first_candidates)
    dict_of_scores_second = create_dict_score_sgs(second_candidates)
    # dictionaries with key = sgRNA seq, value= total score (not including scores below threshold)
    dict_of_scores_first_threshold = create_dict_score_threshold(first_candidates, crispys_threshold)
    dict_of_scores_second_threshold = create_dict_score_threshold(second_candidates, LOW_THRESHOLD_OMEGA)

    # dictionary, where key=subgroup of genes and value = best sgRNAs limited to the permitted number of sgRNAs.
    sgs_per_internal_node_dic = create_dict_of_sgs_per_subgroup_alternative(subgroups_and_sgs_first, 1, number_of_sgs,
                                                                            dict_of_scores_first,
                                                                            crispys_internal_nodes, bottom_up,
                                                                            dict_of_scores_first_threshold)

    update_dict_of_sgs_per_internal_node(sgs_per_internal_node_dic, subgroups_and_sgs_second, subgroups_and_sgs_first,
                                         number_of_sgs, crispys_internal_nodes, dict_of_scores_second, bottom_up,
                                         dict_of_scores_second_threshold)


    if not set(list(sgs_per_internal_node_dic.keys())) <= crispys_internal_nodes:
        raise Exception("Error: join_output(): The insternal nodes do not match!")


    # get the final candidates and the classification of each sgRNA whether it was obtained from the first run or from
    # the second
    final_candidates, dict_class  = find_final_candidates(sgs_per_internal_node_dic, merged_lst_of_candidates,
                                                           sgRNAs_first)


    if not check_if_has_duplcates(sgs_per_internal_node_dic):
        raise Exception("Error: join_output(): There are duplicate sgRNAs!")

    final_pickle_file_path = os.path.join(family_dir, 'res_in_filtered_lst_final.p')
    pickle_path_for_statistics = os.path.join(family_dir, 'dict_sgRNAs_origin.p')
    path_for_internal_nodes_per_sg = os.path.join(family_dir, 'internal_nodes_per_sg.p')

    #sgRNAs_origin_dict = create_origin_dict(filtered_final_list, sgRNAs_first)
    dict_number_sgs_for_all[4] = len(final_candidates)
    pickle.dump(dict_number_sgs_for_all, open(path_for_dic_number_sgs, 'wb'))
    pickle.dump(dict_class, open(pickle_path_for_statistics, 'wb'))
    pickle.dump(final_candidates, open(final_pickle_file_path, 'wb'))
    pickle.dump(sgs_per_internal_node_dic, open(path_for_internal_nodes_per_sg, 'wb'))



    update_lines_of_output(final_candidates, dict_of_lines_first_run, dict_of_lines_second_run)
    new_file_path = os.path.join(family_dir, 'output_filtered_final.csv')


    file = open(new_file_path, 'w')
    final_content = create_final_content(dict_of_lines_first_run, final_candidates)
    file.write(final_content)
    file.close()

    return
#####################################################################################################


def update_lines_of_output(final_candidates, dict_of_lines_first_run, dict_of_lines_second_run):
    add_new_lines(final_candidates, dict_of_lines_first_run, dict_of_lines_second_run)
    set_of_sgs_to_remove = update_set_of_sg_to_remove(final_candidates, dict_of_lines_first_run)
    for sg in set_of_sgs_to_remove:
        dict_of_lines_first_run.pop(sg)
    return



#########################################################################################################

def update_set_of_sg_to_remove(final_candidates, dict_of_lines_first_run):
    set_of_sgs_to_remove = set()
    for sg in dict_of_lines_first_run:
        is_in_lst = False
        for candidate in final_candidates:
            if candidate.seq == sg:
                is_in_lst = True
                break
        if not is_in_lst:
            set_of_sgs_to_remove.add(sg)
    return set_of_sgs_to_remove


#######################################################################################################

def add_new_lines(final_candidates, dict_of_lines_first_run, dict_of_lines_second_run):

    for candidate in final_candidates:

        if not candidate.seq in dict_of_lines_first_run:

            dict_of_lines_first_run[candidate.seq] = dict_of_lines_second_run[candidate.seq]
    return

#######################################################################################################

def create_final_content(dict_of_lines_first_run, filtered_final_list):
    filtered_final_list.sort(key = lambda cand: sum(list(cand.genes_score_dict.values())), reverse= True)
    lst_of_lines = []
    for candidate in filtered_final_list:
        seq = candidate.seq
        lst_of_lines.append(dict_of_lines_first_run[seq])
    joined = '\n'.join(lst_of_lines)
    return joined

########################################################################################################################


def fill_dict_lines(splitted_lines, dict_of_lines):
    p = re.compile("[\S]+")
    prev_line = None
    prev_splitted_line = None
    for line in splitted_lines:
        if line != "":
            splitted_sgRNA_data = re.split(",", line)

            if not prev_line:
                sgRNA = p.findall(splitted_sgRNA_data[1])[0]
                dict_of_lines[sgRNA] = line
            else:
                if len(prev_splitted_line) != len(splitted_sgRNA_data):
                    raise Exception("fill_dict_lines(): lines have different length!!!")
                for i in range(len(splitted_sgRNA_data)):
                    if len(splitted_sgRNA_data[i]) == 0:
                        splitted_sgRNA_data[i] = prev_splitted_line[i]
                sgRNA = p.findall(splitted_sgRNA_data[1])[0]
                ####
                if sgRNA in dict_of_lines:
                    dict_of_lines[sgRNA] += "\n"+ (",".join(splitted_sgRNA_data))
                ####
                else:
                    dict_of_lines[sgRNA] = ",".join(splitted_sgRNA_data)

            prev_line = ",".join(splitted_sgRNA_data)
            prev_splitted_line = splitted_sgRNA_data

            # exists = check_if_exists(lst_candiates, sgRNA)
            # if not exists:
            #     continue
            # dict_of_lines[sgRNA] = line
    return


#######################################################################################################################

def check_if_exists(lst_candidates, sg):
    exists = False
    for candidate in lst_candidates:
        if candidate.seq == sg:
            exists = True
            break
    return exists



##########################################################################################################################
# Input: a list of sgRNAs as candidate instances
# Return Value1: a dictionary of subgroup as a key (exactly the subset of genes to which the sgRNA bounds) and
#   # a list of corresponding candidates as values
# Return Value2: set of sgRNAs (based on the input list)
def create_subgoups(candidates):
    dic = {}
    set_of_sgs = set()
    for candidate in candidates:
        key = tuple(sorted(list(candidate.genes_score_dict.keys())))
        if not key in dic:
            dic[key] = [candidate]
        else:
            dic[key].append(candidate)
        set_of_sgs.add(candidate.seq)
    return dic, set_of_sgs

#########################################################################################################################
# Param1: dictionary of lines from the csv output (the first run)
# Param2: dictionary of lines from the csv output (the second run)
# Return value: For both the runs together- a list of lists, where each list is [sgRNA seq, total score, genes, positions]
# # where position represents a dictionary (key = gene, value = a list of (position, strand (True for '+'/False for '-'))

def create_lst_of_remained(dict_of_lines_first_run, dict_of_lines_second_run):
    list_remained = []
    new_dict_lines = dict_of_lines_second_run.copy()
    new_dict_lines.update(dict_of_lines_first_run)
    p = re.compile("[\S]+")
    for sg in new_dict_lines:
        #####################################################
        splitted_sgRNA_multi_data = re.split("\n", new_dict_lines[sg])
        genes = set()
        genes_and_positions = {}
        for i in range(len(splitted_sgRNA_multi_data)):
            splitted_sgRNA_data = re.split(",", splitted_sgRNA_multi_data[i])
            if i == 0:
                sgRNA = p.findall(splitted_sgRNA_data[1])[0]
                score = p.findall(splitted_sgRNA_data[2])[0]
            gene = p.findall(splitted_sgRNA_data[3])[0]
            genes.add(gene)
            position = p.findall(splitted_sgRNA_data[7])[0]
            if position[-1] == '-':
                if not gene in genes_and_positions:
                    genes_and_positions[gene] = [(int(position[:-1]), False)]
                else:
                    genes_and_positions[gene].append((int(position[:-1]), False))
            else:
                if not gene in genes_and_positions:
                    genes_and_positions[gene] = [(int(position[:-1]), True)]
                else:
                    genes_and_positions[gene].append((int(position[:-1]), True))
        list_remained.append([sgRNA, score, genes, genes_and_positions])

    return list_remained

#################################################################################################################
# Input: list of sgRNAs elements: [sgRNA, score, genes, {gene: [(position, True/False)}]
# Return value: key = sgRNA seq, value = (genes, {gene: [(position, True/False)}, total score)
def create_pos_dict(list_of_sgRNAs_remained):#, cds_dic, list_of_sg):

    dic = {}
    for sgRNA_entry in list_of_sgRNAs_remained:
        sgRNA =  sgRNA_entry[0]
        score = sgRNA_entry[1]
        genes = sgRNA_entry[2]
        positions = sgRNA_entry[3]

        dic[sgRNA] = (genes, positions, float(score))
    return dic



#######################################################################################################################
# Param1: combined dictionary from both the runs: key = subgroup of genes, value = set of sgRNAs.
# Param2: dictionary of cds: key = gene, value = list of exons
# Param3: merged list of sgRNA candidates (from both the runs)
# Param4: dictionary where: key = sgRNA seq, value = (genes, {gene: [(position, True/False)}, total score)
# Param5: threshold distance (for overlap).
# Param6: Crispys threshold (omega)

def filter_same_positions_for_joined(sets_of_sgRNAs, dict_of_cds, proccessed_sg_list, dict_of_pos, thre_distance, crispys_threshold):
    removed_sgs = set()
    threshold_distance = float(thre_distance)
    lst_of_sg_copy = proccessed_sg_list.copy()
    dict_of_indices = create_dict_from_sg_list(lst_of_sg_copy)
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
                        candidate1 = lst_of_sg_copy[dict_of_indices[sgrna1]]
                        candidate2 = lst_of_sg_copy[dict_of_indices[sgrna2]]
                        genes_and_scores_sg1 = candidate1.genes_score_dict
                        genes_and_scores_sg2 = candidate2.genes_score_dict
                        worse_scored = True
                        low_threshold_sgRNA1 = isLowThreshold(genes_and_scores_sg1, crispys_threshold)
                        low_threshold_sgRNA2 = isLowThreshold(genes_and_scores_sg2, crispys_threshold)
                        # if the second sgRNA2 has at least one sgRNA above threshold, and sgRNA1 has none, sgRNA2 should
                        # remain
                        if (low_threshold_sgRNA1) and (not low_threshold_sgRNA2):
                            continue
                        # if we know that sgRNA2 has only score below threshold but sgRNA1 has some above, we would
                        # assign it as worse scored.
                        elif (not low_threshold_sgRNA1) and (low_threshold_sgRNA2):
                            worse_scored = True
                        else:
                            # both are either above threshold or below it.
                            if low_threshold_sgRNA1 and low_threshold_sgRNA2:
                                current_threshold = LOW_THRESHOLD_OMEGA
                            else:
                                current_threshold = crispys_threshold

                            for gene in genes_and_scores_sg2:
                                if genes_and_scores_sg2[gene] > genes_and_scores_sg1[gene]:
                                    if genes_and_scores_sg2[gene] < current_threshold:
                                        continue

                                    worse_scored = False
                                    break

                        overlapping = check_if_overlapping_positions(sgrna1, sgrna2, dict_of_pos, dict_of_cds,
                                                                     threshold_distance, LOW_THRESHOLD_OMEGA, genes_and_scores_sg2)
                        if worse_scored and overlapping:    #It is obvious that if there is such sgRNA, that means that
                                                            # it is from the second run.


                            proccessed_sg_list.remove(candidate2)
                            removed_sgs.add(candidate2.seq)
########################################################################################################################
def isLowThreshold(genes_and_scores_sg, threshold):
    low_threshold = True
    for gene in genes_and_scores_sg:
        if genes_and_scores_sg[gene] >= threshold:
            low_threshold = False
            break
    return low_threshold
########################################################################################################################
# Param1: a dictionary of subgroup as a key (exactly the subset of genes to which the
# # # # # sgRNA bounds) and a list of corresponding candidates as values. This is for the first run
# Param2: the same as param1 but for the second run
# Param3: set of sgRNAs from the first run
# Return value: a combined dictionary of subgroups and the corresponding sgRNAs (with no duplicates)

def merge_subgroups(subgroups_and_sgs_first, subgroups_and_sgs_second, first_sgRNAs):
    merged_final_dic = {}
    for subgroup in subgroups_and_sgs_first:
        merged_final_dic[subgroup] = set([candidate.seq for candidate in subgroups_and_sgs_first[subgroup]])
        if subgroup in subgroups_and_sgs_second:
            candidates_of_second_sub = subgroups_and_sgs_second[subgroup]
            merged_final_dic[subgroup] = \
                merged_final_dic[subgroup].union(
                    set([candidate.seq for candidate in candidates_of_second_sub if not candidate.seq in first_sgRNAs]))
    for subgroup in subgroups_and_sgs_second:
        if not subgroup in subgroups_and_sgs_first:
            merged_final_dic[subgroup] = set([candidate.seq for candidate in subgroups_and_sgs_second[subgroup] if not \
                candidate.seq in first_sgRNAs])

    return merged_final_dic

###############################################################################################################################
# Param1: list of candidate instances from the first run
# Param2:  list of candidate instances from the second run
# Param3: set of sgRNAs from the first run
# Return value: combined list of sgRNAs (no duplicates)

def merge_list_of_candidates(first_candidates, second_candidates, first_sgs):
    lst_final = first_candidates.copy()
    for candidate in second_candidates:
        if candidate.seq in first_sgs:
            continue
        lst_final.append(candidate)
    return lst_final


########################################################################################################################

def update_candidates_list(sgRNAs_set, merged_lst_of_candidates, first_sgRNAs):
    new_list_of_candidates = []
    for candidate in merged_lst_of_candidates:
        if candidate.seq in sgRNAs_set:
            if not first_sgRNAs:
                new_list_of_candidates.append(candidate)
            else:
                if not candidate.seq in first_sgRNAs:
                    new_list_of_candidates.append(candidate)
    return new_list_of_candidates

#######################################################################################################################
# Param1: dictionary for the first run, where key=subgroup of genes and value = best sgRNAs limited to the
# permitted number of sgRNAs.
# Param2: a dictionary of subgroup as a key (exactly the subset of genes to which the
# sgRNA bounds) and a list of corresponding candidates as values (second run)
# Param3:  a dictionary of subgroup as a key (exactly the subset of genes to which the
# sgRNA bounds) and a list of corresponding candidates as values (first run)
# Param4: number of sgRNAs per internal node
# Param5: set of internal nodes represented as tuples of genes in a sorted order.
# Param6: a dictionary with key = sgRNA seq, value= total score
# Param7: 1 if bottom-up. 0 otherwise.
# Param8: dictionary for the second run- key = sgRNA seq, value= total score (not including scores below threshold)
# This dictionary is for the second run.
# the function updates the internal nodes with new sgRNAs (from the second run).
def update_dict_of_sgs_per_internal_node(sgs_per_internal_node_dic, subgroups_and_sgs_second,
                                         subgroups_and_sgs_first, number_of_sgs, internal_nodes, dict_of_scores,
                                         bottom_up_indicator, dict_of_scores_second_threshold):
    set_of_used_sgrnas = set()
    set_of_all_possible_subgroups = set(subgroups_and_sgs_first.keys()).union(set(subgroups_and_sgs_second.keys()))
    set_of_all_possible_subgroups = set_of_all_possible_subgroups.union(internal_nodes)
    dict_of_key_subsets = create_dict_with_dummy_values(set_of_all_possible_subgroups)
    print("dict_of_key_subsets: ", dict_of_key_subsets)
    print("**")
    dict_of_subsets, sorted_sets = \
        create_dict_of_subsets(dict_of_key_subsets, bottom_up = bottom_up_indicator, top_down = not bottom_up_indicator)
    print("dictionary of sorted: ", dict_of_subsets)
    for gene_set in sorted_sets:
        if not gene_set in internal_nodes:
            continue
        if gene_set in sgs_per_internal_node_dic:
            num_of_sgs_from_first_run = len(sgs_per_internal_node_dic[gene_set])
        else:
            num_of_sgs_from_first_run = 0
        num_to_add = number_of_sgs - num_of_sgs_from_first_run
        print("*****\n")
        print("Main set: ", gene_set)

        if num_to_add > 0:

            sgrnas = find_sgrnas_to_add(dict_of_subsets[gene_set], subgroups_and_sgs_second, dict_of_scores,
                                        set_of_used_sgrnas, dict_of_scores_second_threshold)


            for i in range(min(num_to_add, len(sgrnas))):
                if not gene_set in sgs_per_internal_node_dic:
                    sgs_per_internal_node_dic[gene_set] = {sgrnas[i]}
                    set_of_used_sgrnas.add(sgrnas[i])
                else:
                    sgs_per_internal_node_dic[gene_set].add(sgrnas[i])
                    set_of_used_sgrnas.add(sgrnas[i])
    return

##################################################################################################################
# Param1: subsets of the current set of genes
# Param2: a dictionary of subgroup as a key (exactly the subset of genes to which the
# # sgRNA bounds) and a list of corresponding candidates as values (second run).
# Param3: a dictionary with key = sgRNA seq, value= total score
# Param4: set of used sgRNAs (they have been already added)
# Param5: dictionary for the second run- key = sgRNA seq, value= total score (not including scores below threshold)
# second run.
# Return Value: a list of sgRNAs which correspond to the subsets of genes of a particular set of genes. Sorted according to the score
# in a descending order (the best are the first)
def find_sgrnas_to_add(subsets, subgroups_candidates, dict_of_scores, set_of_used, dict_of_scores_second_threshold):
    sgs = []
    for subgroup in subsets:
        sorted_gene_set = tuple(sorted(list(subgroup)))
        if not sorted_gene_set in subgroups_candidates:
            continue
        print("subset is: ", sorted_gene_set)
        sgs.extend([candidate.seq for candidate in subgroups_candidates[sorted_gene_set] if not candidate.seq in set_of_used])
    #sgs.sort(reverse= True, key = lambda x : dict_of_scores_second_threshold[x])
    print("sgRNAs: ", sgs)
    sgs = sort_according_to_score(sgs, dict_of_scores, dict_of_scores_second_threshold)
    print("sgRNAs after soring", sgs)
    return sgs
#######################################################################################################################


def find_final_candidates(sgs_per_internal_node_dic, merged_lst_of_candidates, sgRNAs_first):
    final_list = []
    dict_class = {}
    set_of_sgrnas = set()
    for gene_set in sgs_per_internal_node_dic:
        set_of_sgrnas = set_of_sgrnas.union(sgs_per_internal_node_dic[gene_set])
    for candidate in merged_lst_of_candidates:
        if candidate.seq in set_of_sgrnas:
            if candidate.seq in sgRNAs_first:
                dict_class[candidate.seq] = 1
            else:
                dict_class[candidate.seq] = 2
            final_list.append(candidate)
    return final_list, dict_class



#######################################################################################################################

def check_if_has_duplcates(sgs_per_internal_node_dic):
    set_of_sg = set()
    counter = 0
    for gene_set in sgs_per_internal_node_dic:
        counter += len(sgs_per_internal_node_dic[gene_set])
        set_of_sg = set_of_sg.union(sgs_per_internal_node_dic[gene_set])

    return len(set_of_sg) == counter

########################################################################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filters off")###################################################################
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--second_dir", "-s", help = "directory of second run")
    parser.add_argument("--family", "-f", help = "family")
    parser.add_argument("--number_of_sgs", "-n", help = "number of sgRNAs per internal node")
    parser.add_argument("--thr_distance", "-y", help = "threshold distance")
    parser.add_argument("--bottom_up", "-b", help = "specify if it is bottom up")
    parser.add_argument("--crispys_threshold", "-c", help = 'crispys threshold')


    args = parser.parse_args()
    output = args.output
    second_run_dir = args.second_dir
    family = args.family
    number_of_sgs = args.number_of_sgs
    threshold_distance = args.thr_distance
    bottom_up = args.bottom_up
    crispys_threshold = args.crispys_threshold

    join_output(output, second_run_dir, family, int(number_of_sgs), float(threshold_distance),
                int(bottom_up), float(crispys_threshold))