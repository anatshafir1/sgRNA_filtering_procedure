import argparse
import numpy as np
from create_dictionaries_and_pickle import *
from subprocess import Popen
from global_variables import *


#Main finction, which excludes sgRNAs, which are off- target
#Param1: output directory (Directory of all the families
#Return: No return value. A direcory with output files will be created

def main(target_dir, threshold, threshold_distance, number_of_sgs, bottom_up, method, crispys_threshold, relative):   #target_dir must be the directory of all the families
    target_dir_annotation_file = os.path.join(MAIN_FILES_DIR, "annotation")
    genome_data_path = os.path.join(MAIN_FILES_DIR, "ath.con")
    dic_genome = create_chromosome_seq_file(genome_data_path)
    family_dic_path = os.path.join(MAIN_FILES_DIR, ORGANISM_SHORT_NAME +"_genefamily_dict_without_mt_cp_genes.p")
    annotation_file_source = os.path.join(MAIN_FILES_DIR, "annotation.ath.csv")
    pickle_genome_file_path = os.path.join(MAIN_FILES_DIR, "genome_dic.p")
    if not os.path.exists(pickle_genome_file_path):
        pickle.dump(dic_genome, open(pickle_genome_file_path, "wb"))
        print("Created the chromosome pickle file ....")
    list_of_dirs = os.listdir(target_dir)

    if not os.path.exists(target_dir_annotation_file):
        os.makedirs(target_dir_annotation_file)

    pickle_annotation_file_path = os.path.join(target_dir_annotation_file, "annotation.p")
    if not os.path.exists(pickle_annotation_file_path):
        parse_annotation_file(annotation_file_source, pickle_annotation_file_path)
        print("Parsed the annotation file ....")

    family_genes_dic = load_pickle_files(family_dic_path)
    annotation_data = load_pickle_files(pickle_annotation_file_path)

    path_for_exons_dic = os.path.join(MAIN_FILES_DIR, 'dict_exons.p')
    dic_exons = getExonsDictionary(annotation_data)
    pickle.dump(dic_exons, open(path_for_exons_dic, "wb"))
    print("Written the dic_exons.p file ....")
    print("Starting the iteration over families!")

    for dir in list_of_dirs:
        full_dir_path = os.path.join(target_dir, dir) # source directory
        if not os.path.isdir(full_dir_path):    # iterate only over directories (of families)
            continue
        family_genes = family_genes_dic[dir]
        annotation_path_for_family = os.path.join(full_dir_path, 'annotation.p')
        # create annotation data structure for each family
        annotation_dic_for_family = {}
        for gene in family_genes:
            annotation_dic_for_family[gene] = annotation_data[gene]
        pickle.dump(annotation_dic_for_family, open(annotation_path_for_family, "wb"))


        job = create_job(full_dir_path, dir, annotation_path_for_family,
                         pickle_genome_file_path, full_dir_path, threshold, threshold_distance,
                         number_of_sgs, bottom_up, method, crispys_threshold, relative)
        Popen(["qsub", "-p", "-1", job])
    return

#######################################################################################################################
def create_job(full_dir_path, dir, annotation_file_path, dic_genome_path, target_dir, threshold, threshold_distance,
               number_of_sgs, bottom_up, method, crispys_threshold, relative):
    python_file = os.path.join(SCRIPSTS_DIR, "filter_sgRNAs.py")
    second_run = 0
    if full_dir_path != target_dir:
        second_run = 1

    with open(full_dir_path + "/job_python" + str(second_run)+ ".sh", "w") as handle:
        handle.write(createHeaderJob(full_dir_path, "sg_Filtering" + str(second_run)))
        handle.write(PYTHON_MODULE)
        handle.write("python "+python_file+" -o " + full_dir_path +" -d "+dir+" -a "+\
                     annotation_file_path+ " -g "+ dic_genome_path+ " -t "+target_dir+ " -w "+threshold +" -y "+
                     threshold_distance+" -n "+number_of_sgs+ " -b "+bottom_up+ " -m "+method+ " -c "+
                     crispys_threshold+
                     " -r " + relative +"\n")

    return full_dir_path + "/job_python"+ str(second_run)+".sh"


#######################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filters off")#########################################################
    parser.add_argument("--output", "-o", help = "output directory") # full path
    parser.add_argument("--threshold", "-w", help = "threshold for off- target")
    parser.add_argument("--threshold_distance", "-y", help = "threshold for overlapping")
    parser.add_argument("--number_of_sgs", "-n", help = "number of sgRNAs per internal node")
    parser.add_argument("--bottom_up", "-b", help = "specify if the filtering and addition is bottom up or top down")
    parser.add_argument("--method", "-m", help = "if running just once- 0. Otherwise- 1")   #need to remove this option !!!
    parser.add_argument("--crispys_threshold", "-c", help = "crispys threshold")
    parser.add_argument("--relative", "-r", help = 'relative')

    args = parser.parse_args()
    output = args.output
    threshold = args.threshold
    threshold_dis = args.threshold_distance
    number_of_sgs = args.number_of_sgs
    bottom_up = args.bottom_up  # 0 or 1
    method = args.method
    crispys_threshold = args.crispys_threshold
    relative = args.relative

    main(output, threshold, threshold_dis, number_of_sgs, bottom_up, method, crispys_threshold, relative)


