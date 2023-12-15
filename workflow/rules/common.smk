import glob
import os


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


SIMULATION_OUTPUT = expand("results/simulate_fastas/sample_{SID}_{READS}.fasta", 
                            SID = sample_names,
                            READS = READS)



### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Directory containing index; used in case of certain aligners
INDEX_DIR = config["indices"]

# Get input fastq files for first step
if PE:

    FASTQ_FILES = [config["read_1"], config["read_2"]]

else:

    FASTQ_FILES = list(config["read_1"])


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

    return target