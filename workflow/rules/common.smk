### Simulation output

if config["read_2"]:

    READS = ['1', '2']

    SIMULATION_PARAMS = config["simulate_fastas_params"]

else:

    READS = ['1']

    SIMULATION_PARAMS = config["simulate_fastas_params"] + " --singleend"

def generate_formatted_list(n):
    formatted_list = ["{:02d}".format(i) for i in range(1, n + 1)]
    return formatted_list

sample_names = generate_formatted_list(config["number_of_replicates"])

SIMULATION_OUTPUT = expand("results/simulate_fastas/sample_{S_ID}_{READS}", 
                            SID = sample_names,
                            READS = READS)