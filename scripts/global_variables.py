import pickle
SCRIPSTS_DIR = "/path/to/scripts/directory"
MAIN_FILES_DIR = "/path/to/directory/of/that/contains/families"
PYTHON_MODULE = "module load your/python/version\n"
DIR_NAME_BWA = "BWA"
LOW_THRESHOLD_OMEGA = 0.45
THRESHOLD_NUM_SAME_POSITONS_PER_GENE = 0.9
THRESHOLD_FRAC_WORSE_SCORE = 0.5
SCORE_FUNCTION = "cfd_eps_dict.p"
BWA_MODULE = "module load /your/bwa/version\n" # updated according to Power
ORGANISM_SHORT_NAME = "arabidopsis"
                                            # (the previous command was relevant for Jekyl and lecs)
########################################################################

def createHeaderJob(path, job_name, ncpu = 1):
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q your_queue\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N "+ job_name+"\n"
    text += "#PBS -e " + path + "/"+job_name+".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name +".OU"+ "\n"
    text += "#PBS -l select=ncpus="+ str(ncpu)+ ":mem=2gb\n"
    return text
#######################################################################33
#A function which loads a Pickle file.
#Input: pickle file path.
#Output: dictionary

def load_pickle_files(pickle_file):
    with open(pickle_file, 'rb') as handle:
        dic = pickle.load(handle)
    handle.close()
    return dic
