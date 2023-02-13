from filter_sgRNAs import *
from  sgRNA_classification import *
from subprocess import Popen

# This function should join the results of two Crispys runs
# Param1: the directory of all families
# Param2: the number of sgRNAs per internal node
# Param3: 1 if bottom-up. 0 otherwise.
# Param4: Crispys threshold #TODO: which one?
# No return value: runs for each family the joined filtering procedure

def combine_results(output_dir, number_of_sgs, threshold_distance, bottom_up, crispys_threshold):
    lst_of_dirs = os.listdir(output_dir)
    for dir in lst_of_dirs:
        family_dir = os.path.join(output_dir, dir)
        if not os.path.isdir(family_dir):
            continue
        path_dir_second_run = os.path.join(family_dir, 'second_run')
        if not os.path.exists(path_dir_second_run):
            continue

        job = create_job_for_join(family_dir, path_dir_second_run, dir, str(number_of_sgs),
                                  str(threshold_distance), str(bottom_up), crispys_threshold)
        Popen(["qsub", "-p", "-1", job])

########################################################################################################################
def create_job_for_join(family_dir, path_dir_second_run, dir, number_of_sgs,
                        threshold_distance, bottom_up, crispys_threshold):

    python_file = os.path.join(SCRIPSTS_DIR, "join_outputs.py")
    with open(family_dir + "/job_python_join.sh", "w") as handle:
        handle.write(createHeaderJob(family_dir, "join_two_runs", ncpu=1))
        handle.write("cd " + family_dir + "\n")
        handle.write(PYTHON_MODULE)
        handle.write("python "+python_file+" -o " + family_dir +" -s "+path_dir_second_run+" -f "+\
                     dir+ " -y "+threshold_distance+" -n "+number_of_sgs+ " -b "+bottom_up+" -c "\
                     +crispys_threshold+"\n")
    return family_dir + "/job_python_join.sh"
########################################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filters off")
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--number_of_sgs", "-n", help = "number of sgRNAs per internal node")
    parser.add_argument("--thr_distance", "-t", help = "threshold_distance")
    parser.add_argument("--bottom_up", "-b", help = "threshold_distance")
    parser.add_argument("--crispys_threshold", "-c", help = "crispys threshold")

    args = parser.parse_args()
    output = args.output
    number_of_sgs = args.number_of_sgs
    threshold_distance = args.thr_distance
    bottom_up = args.bottom_up
    crispys_threshold = args.crispys_threshold

    print("running")

    combine_results(output, int(number_of_sgs), float(threshold_distance), int(bottom_up), crispys_threshold)
