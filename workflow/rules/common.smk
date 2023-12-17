import glob
import os
import math


### Simulation output

if config["read_2"]:

    READS = ['1', '2']

    SIMULATION_PARAMS = config["simulate_fastas_params"]

    PE = True

else:

    READS = ['1']

    SIMULATION_PARAMS = config["simulate_fastas_params"] + " --singleend"

    PE = False


def generate_formatted_list(n):
    formatted_list = ["{:02d}".format(i) for i in range(1, n + 1)]
    return formatted_list


sample_names = generate_formatted_list(config["number_of_replicates"])

if PE:

    SIMULATION_OUTPUT = expand("results/simulate_fastas/sample_{SID}_{READS}.fasta", 
                                SID = sample_names,
                                READS = READS)

else:

    SIMULATION_OUTPUT = expand("results/simulate_fastas/sample_{SID}.fasta", 
                                SID = sample_names)


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Make life easier for users and catch if they add a '/' at the end of their path
# to alignment indices. If so, remove it to avoid double '/' 
if config["indices"].endswith('/'):
    INDEX_PATH = str(config["indices"])
else:
    INDEX_PATH = str(config["indices"]) + "/"

# Get input fastq files for first step
if PE:

    FASTQ_FILES = [config["read_1"], config["read_2"]]

else:

    FASTQ_FILES = [config["read_1"]]


# Check if fastq files are gzipped
is_gz = False


test_gz = any(path.endswith('.fastq.gz') for path in FASTQ_FILES)
is_gz = any([is_gz, test_gz])




### SALMON HELPERS

# Libtype string
if PE:

    if config["strandedness"] == "reverse":

        LIBTYPE = config["directionality"] + "SR"
    
    if config["strandedness"] == "yes":

        LIBTYPE = config["directionality"] + "SF"
    
    if config["strandedness"] == "no":

        LIBTYPE = config["directionality"] + "U"

else:

    if config["strandedness"] == "reverse":

        LIBTYPE = "SR"
    
    if config["strandedness"] == "yes":

        LIBTYPE = "SF"
    
    if config["strandedness"] == "no":

        LIBTYPE = "U"


### Target rule

def get_target_output():

    target = SIMULATION_OUTPUT

    if config["run_fastqc"]:

        target.append(expand("results/fastqc/read{read}.html", read = READ_NAMES))

    if PE:

        target.append(expand("results/convert_to_fastq/sample_{sample}_{read}.fastq.gz", sample = sample_names, read = READ_NAMES))

    else:

        target.append(expand("results/convert_to_fastq/sample_{sample}.fastq.gz", sample = sample_names))

    return target


### FASTA file splitting

# Number of temporary fasta files to make
NUMBER_SPLIT = math.ceil(config["library_size"]/32)

# Number of digits in number of temp files; for inferring file name
num_digits = len(str(NUMBER_SPLIT))

# IDS in split FASTA file sample names
SPLIT_IDS = [str(i).zfill(num_digits) for i in range(0, NUMBER_SPLIT)]