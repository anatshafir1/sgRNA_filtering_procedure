import crista_search_offtargets
from filter_sgRNAs import *
from global_variables import*

# execte bwa for a given family
# Param1: a full path for the family directory
# Param2: Cripys threshold for the current run
# Param3: 1 if it is the second run (with a lower threshold). Otherwise, 0.
# Return value: No return.
def execute_bwa(full_dir_path, crispys_threshold, second_run):
    if second_run:
        family = os.path.basename(os.path.normpath(full_dir_path))
        lst_of_all_sg_pickle = os.path.join(full_dir_path, "second_run", family, 'res_in_lst.p')
    else:
        lst_of_all_sg_pickle = os.path.join(full_dir_path, 'res_in_lst.p')
    with open(lst_of_all_sg_pickle, 'rb') as handle:
        lst_sgRNAs = pickle.load(handle)  #A list of lists of sgRNAs (divided into sets)
    handle.close()
    if lst_sgRNAs == []:
        return

    proccessed_sg_list = concatenate_list(lst_sgRNAs)    # concatenates all the subsets
    delete_singletons(proccessed_sg_list, crispys_threshold)
    leave_only_relevant_sgRNA(proccessed_sg_list)    #first filtering of sgRNAs
    # create a directory where all the bwa results will be stored
    # for each sgRNA create a directory in the following format sgRNA_{sgRNA}
    file_suffixes_to_del = ["sgrna.fa", "sgrna.fa.out.sai", "sgrna.fa.out.sam", "sgrna.fa.out.table"]
    family_new_dir_path_sgRNA = os.path.join(full_dir_path, DIR_NAME_BWA, "off_target_results")
    if not os.path.exists(family_new_dir_path_sgRNA):
        os.makedirs(family_new_dir_path_sgRNA)
    for candidate in proccessed_sg_list:
        sgRNA = candidate.seq
        #family_new_dir_path_sgRNA = os.path.join(full_dir_path, DIR_NAME_BWA, "sgRNA_" + sgRNA)
        #if not os.path.exists(family_new_dir_path_sgRNA):
            #os.makedirs(family_new_dir_path_sgRNA)
        if find_crista_output_filename(family_new_dir_path_sgRNA, sgRNA):
            continue
        crista_search_offtargets.main(sgRNA, family_new_dir_path_sgRNA)
        for file_to_del in file_suffixes_to_del:
            fileName = os.path.join(full_dir_path, DIR_NAME_BWA, "off_target_results" + file_to_del)
            if os.path.exists(fileName):
                os.remove(fileName)
            ## due to some wierd behaviour of the script after version updates

    print("Job ended successfully!") # useful for testing
    return

########################################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filters off")
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--threshold", "-t", help = "threshold")
    parser.add_argument("--second_run", "-s", type = int, help = "is it a second run?")

    args = parser.parse_args()
    output = args.output
    threshold = args.threshold
    second_run = args.second_run

    execute_bwa(output, float(threshold), second_run)