from subprocess import Popen
import os
from global_variables import *
import argparse

def run_Crispys(output_dir, second_run, Crispys_threshold):
    famlies = os.listdir(output_dir)
    for family in famlies:
        target_dir = None
        full_family_path = os.path.join(output_dir, family)
        if not os.path.isdir(full_family_path):
            continue
        if not second_run:
            target_dir = full_family_path
        else:
            second_run_dir = os.path.join(full_family_path, "second_run", family)
            if not os.path.exists(second_run_dir):
                continue
            target_dir = second_run_dir
            if LOW_THRESHOLD_OMEGA != Crispys_threshold:
                raise Exception("ERROR!! run_Crispys(): Crispys threshold of second run differs from the global LOW_THRESHOLD_OMEGA!")
        if not target_dir:
            raise Exception("ERROR!! run_Crispys(): target_dir was not assigned to any value!!")
        createCrispysJob(target_dir, Crispys_threshold, second_run, family)
#######################################################################################################################

def createCrispysJob(target_dir, Crispys_threshold, second_run, family):
    #python_script_path = "/groups/itay_mayrose/galhyams/Crispys_21_02/call_MULTICRISPR_Wrapper.py"
    python_script_path = "/groups/itay_mayrose/udiland/CRISPys-master/Stage0.py"
    fasta_file = os.path.join(target_dir, family +".txt")
    header = createHeaderJob(target_dir, "Crispys_"+ str(second_run))
    command = "module load mafft/mafft7149\n"
    command += "module load python/anaconda3-5.0.0\n"
    command += "cd "+ target_dir +"\n"

    program_command = "python "+ python_script_path +" "+ fasta_file + " "+ target_dir + " --alg E --where_in_gene 0.67 --t 1 --v "\
                      + str(Crispys_threshold)+" --i 200 --s cfd_funct\n"
    job_content = header+ command + program_command
    job_path = os.path.join(target_dir, "Crispys"+ str(second_run)+".sh")
    with open (job_path, 'w') as handle:
        handle.write(job_content)
    Popen(["qsub", "-p", "-1", job_path])

########################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="filters off")
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--second_run", "-s", type = int, help="second run of Cripys")
    parser.add_argument("--crispys_threshold", "-c", help = "crispys_threshold")


    args = parser.parse_args()
    output = args.output
    second_run = args.second_run
    crispys_threshold = args.crispys_threshold

    run_Crispys(output, second_run, crispys_threshold)


