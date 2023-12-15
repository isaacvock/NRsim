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


SIMULATION_OUTPUT = expand("results/simulate_fastas/sample_{S_ID}_{READS}", 
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
fastq_paths = config["samples"]

is_gz = False

for p in fastq_paths.values():

    fastqs = sorted(glob.glob(f"{p}/*.fastq*"))
    test_gz = any(path.endswith('.fastq.gz') for path in fastqs)
    is_gz = any([is_gz, test_gz])




### SALMON HELPERS

# Trimmed fastq file paths, used as input for aligners

def get_fastq_r1(wildcards):

    return expand("results/trimmed/read.1.fastq", SID = wildcards.sample)

def get_fastq_r2(wildcards):

    if PE:

        return expand("results/trimmed/read.2.fastq", SID = wildcards.sample)

    else:

        return ""


# Libtype string if using salmon
if config["PE"]:

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


