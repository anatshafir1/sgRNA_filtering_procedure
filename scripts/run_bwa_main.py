from subprocess import Popen
import os
import argparse
from global_variables import*

# This script runs BWA for each gene family
def run_main(families_dir, threshold, second_run = 0):
    families = os.listdir(families_dir)
    for family in families:
        full_dir = os.path.join(families_dir, family)
        if not os.path.isdir(full_dir):
            continue
        job = create_bwa_job(full_dir, threshold, second_run)
        Popen(["qsub", "-p", "-1", job])
    return

def create_bwa_job(full_dir_path, threshold, second_run):
    python_file = os.path.join(SCRIPSTS_DIR, "run_bwa.py")
    with open(full_dir_path + "/bwa_job" + str(second_run)+".sh", "w") as handle:
        handle.write(createHeaderJob(full_dir_path, "BWA"+ str(second_run)))
        handle.write("cd " + full_dir_path + "\n")
        handle.write(PYTHON_MODULE)
        handle.write(BWA_MODULE)
        handle.write("python "+python_file+" -o " + full_dir_path +" -t "+ threshold+ " -s " + str(second_run)+"\n")
    return full_dir_path + "/bwa_job"+ str(second_run)+".sh"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filters off")
    parser.add_argument("--output", "-o", help = "output directory") # full path
    parser.add_argument("--threshold", "-t", help="threshold")
    parser.add_argument("--second_run", "-s", type = int, help = "second run")

    args = parser.parse_args()
    output = args.output
    threshold = args.threshold
    second_run = args.second_run
    run_main(output, threshold, second_run)
