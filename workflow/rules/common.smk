import glob
import os
import math


### Simulation output

SIM_NAMES = list(config['simulation_parameters']['number_of_replicates'].keys())


if config["read_2"]:

    READS = ['1', '2']

    PE_input = True
  

    SIMULATION_PARAMS = config["simulate_fastas_params"]

else:

    READS = ['1']

    PE_input = False

    SIMULATION_PARAMS = config["simulate_fastas_params"]


PE = config["PE"]


def generate_formatted_list(n):
    formatted_list = ["{:02d}".format(i) for i in range(1, n + 1)]
    return formatted_list


sample_names = generate_formatted_list(config["number_of_replicates"])

if config["simulation_parameters"]:

    def get_simulation_output(wildcards):

        if config["simulation_parameters"]["pe"][str(wildcards.sim)]:

            READS_SIM = ['1', '2']


            SIMULATION_DATASET = expand("results/simulate_fastas/{SIM}/sample_{SID}_{READS}.fasta", 
                                        SIM = wildcards.sim,
                                        SID = sample_names,
                                        READS = READS_SIM)

        else:

            SIMULATION_DATASET = expand("results/simulate_fastas/{SIM}/sample_{SID}.fasta", 
                                        SIM = wildcards.sim,
                                        SID = sample_names)

        return SIMULATION_DATASET

else:

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
if PE_input:

    FASTQ_FILES = [config["read_1"], config["read_2"]]

else:

    FASTQ_FILES = [config["read_1"]]


# Check if fastq files are gzipped
is_gz = False


test_gz = any(path.endswith('.fastq.gz') for path in FASTQ_FILES)
is_gz = any([is_gz, test_gz])




### SALMON HELPERS

# Libtype string
if PE_input:

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

        target.append(expand("results/fastqc/read{read}.html", read = READS))


    if config["simulation_parameters"]:

        for s in SIM_NAMES:

            if config["simulation_parameters"]["pe"][s]:

                target.append(expand("results/convert_to_fastq/{sim}/sample_{sample}_{read}.fastq.gz", sim = s, 
                                                                                                    sample = generate_formatted_list(config["simulation_parameters"]["number_of_replicates"][s]), 
                                                                                                    read = ['1', '2']))

            else:

                target.append(expand("results/convert_to_fastq/{sim}/sample_{sample}.fastq.gz", sim = s, 
                                                                                                    sample = generate_formatted_list(config["simulation_parameters"]["number_of_replicates"][s])))


    else:

        if PE:

            target.append(expand("results/convert_to_fastq/sample_{sample}_{read}.fastq.gz", sample = sample_names, read = READS))

        else:

            target.append(expand("results/convert_to_fastq/sample_{sample}.fastq.gz", sample = sample_names))

    return target


### FASTA file splitting

# Number of temporary fasta files to make
    # Hardcoding for now to maximum amount of threads allowed by rule
    # Could get fancier or allow more user specification here, but this is probably sufficient
NUMBER_SPLIT = 32

# Number of digits in number of temp files; for inferring file name
num_digits = len(str(NUMBER_SPLIT))

# IDS in split FASTA file sample names
SPLIT_IDS = [str(i).zfill(num_digits) for i in range(0, NUMBER_SPLIT)]