from filter_sgRNAs_main import *
import argparse
from subprocess import Popen
from global_variables import *



def run_2_round(output_dir, number_of_sgs, threshold_off_target, threshold_dis, bottom_up, method, crispys_threshold, relative):
    lst_of_dirs = os.listdir(output_dir)
    for dir in lst_of_dirs:
        family_dir = os.path.join(output_dir, dir)
        if not os.path.isdir(family_dir):
            continue
        path_dir_second_run = os.path.join(family_dir, 'second_run')
        if not os.path.exists(path_dir_second_run):
            continue
        annotation_path_for_family = os.path.join(family_dir, "annotation.p")
        pickle_genome_file_path = os.path.join(MAIN_FILES_DIR, "genome_dic.p")

        job = create_job(os.path.join(path_dir_second_run, dir), dir,
                         annotation_path_for_family, pickle_genome_file_path,
                         family_dir ,threshold_off_target, threshold_dis, number_of_sgs, bottom_up, method,
                         crispys_threshold, relative)

        Popen(["qsub", "-p", "-1", job])





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filters off")###################################################################
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--threshold", "-w", help = "threshold for determing if a specific sgRNA is off-target")
    parser.add_argument("--threshold_distance", "-y", help = "threshold for overlapping")
    parser.add_argument("--number_of_sgs", "-n", help = "number of sgRNAs per internal node")
    parser.add_argument("--bottom_up", "-b", help = "bottom up or not")
    parser.add_argument("--method", "-m", help = "use the default or the alternative")
    parser.add_argument("--crispys_threshold", "-c", help = "crispys threshold")
    parser.add_argument("--relative", "-r", help = "relative")

    args = parser.parse_args()
    output = args.output
    threshold = args.threshold
    threshold_dist = args.threshold_distance
    number_of_sgs = args.number_of_sgs
    bottom_up = args.bottom_up
    method = args.method
    crispys_threshold = args.crispys_threshold
    relative = args.relative


    run_2_round(output, number_of_sgs, threshold, threshold_dist, bottom_up, method, crispys_threshold, relative)
