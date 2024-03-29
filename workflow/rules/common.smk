import glob
import os
import random
import math


### Simulation output

if config["dataset_specific"]:
    SIM_NAMES = list(config["simulation_parameters"]["library_size"].keys())


if config["read_2"]:
    READS_input = ["1", "2"]

    PE_input = True

    SIMULATION_PARAMS = config["simulate_fastas_params"]

else:
    READS_input = ["1"]

    PE_input = False

    SIMULATION_PARAMS = config["simulate_fastas_params"]


PE = config["PE"]

if PE:
    READS = ["1", "2"]

else:
    READS = ["1"]


def generate_formatted_list(n):
    formatted_list = ["{:02d}".format(i) for i in range(1, n + 1)]
    return formatted_list


sample_names = generate_formatted_list(config["number_of_replicates"])


if PE:
    SIMULATION_OUTPUT = expand(
        "results/simulate_fastas/sample_{SID}_{READS}.fasta",
        SID=sample_names,
        READS=READS,
    )

else:
    SIMULATION_OUTPUT = expand(
        "results/simulate_fastas/sample_{SID}.fasta", SID=sample_names
    )

### Randomly create seed that will be set to ensure reads in pair have same newness status
seeds = {}

if config["dataset_specific"]:

    for sim in SIM_NAMES:
        seeds[sim] = {}

        for samp in sample_names:
            seeds[sim][samp] = 42 # round(random.random() * 10000000) + 1

else:

    for samp in sample_names:
        seeds[samp] = 42 # round(random.random() * 10000000) + 1


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Make life easier for users and catch if they add a '/' at the end of their path
# to alignment indices. If so, remove it to avoid double '/'
if config["indices"].endswith("/"):
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


test_gz = any(path.endswith(".fastq.gz") for path in FASTQ_FILES)
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
    target = []

    if config["run_fastqc"]:
        target.append(expand("results/fastqc/read{read}.html", read=READS))

    if config["dataset_specific"]:
        if PE:
            target.append(
                expand(
                    "results/convert_to_fastq/{sim}/sample_{sample}_{read}.fastq.gz",
                    sim=SIM_NAMES,
                    sample=sample_names,
                    read=READS,
                )
            )
            target.append(
                expand(
                    "results/shuffle_fastq/{sim}/sample_{sample}/read_{read}.fastq.gz",
                    sim=SIM_NAMES,
                    sample=sample_names,
                    read=READS,
                )
            )

        else:
            target.append(
                expand(
                    "results/convert_to_fastq/{sim}/sample_{sample}.fastq.gz",
                    sim=SIM_NAMES,
                    sample=sample_names,
                )
            )
            target.append(
                expand(
                    "results/shuffle_fastq/{sim}/sample_{sample}/read_1.fastq.gz",
                    sim=SIM_NAMES,
                    sample=sample_names,
                )
            )

    else:
        if PE:
            target.append(
                expand(
                    "results/convert_to_fastq/sample_{sample}_{read}.fastq.gz",
                    sample=sample_names,
                    read=READS,
                )
            )
            target.append(
                expand(
                    "results/shuffle_fastq/sample_{sample}/read_{read}.fastq.gz",
                    sample=sample_names,
                    read=READS,
                )
            )

        else:
            target.append(
                expand(
                    "results/convert_to_fastq/sample_{sample}.fastq.gz",
                    sample=sample_names,
                )
            )
            target.append(
                expand(
                    "results/shuffle_fastq/sample_{sample}/read_1.fastq.gz",
                    sample=sample_names,
                )
            )

    return target


### FASTA file splitting

# Number of temporary fasta files to make
# Hardcoding for now to maximum amount of threads allowed by rule
# Could get fancier or allow more user specification here, but this is probably sufficient
NUMBER_SPLIT = 32

# Number of digits in number of temp files; for inferring file name
num_digits = len(str(NUMBER_SPLIT))

# IDS in split FASTA file sample names
SPLIT_IDS = [str(i) for i in range(1, NUMBER_SPLIT + 1)]
## Old SPLIT_IDS for pyfasta
# SPLIT_IDS = [str(i).zfill(num_digits) for i in range(0, NUMBER_SPLIT)]
