rule filter_annotation:
    input:
        gtf="results/modify_annotation/modified_annotation.gtf",
        quantification="results/quant/quant.sf"
    output:
        filtered="results/filter_annotation/filtered_annotation.gtf",
        counts="results/filter_annotation/transcript_read_counts.csv",
    params:
        rscript=workflow.source_path("../scripts/filter_annotation.R"),
        extra=config["filter_annotation_params"]
    conda:
        "../envs/simulate.yaml"
    log:
        "logs/filter_annotation/filter_annotation.log"
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -g {input.gtf} -q {input.quantification} \
        -c {output.counts} -o {output.filtered} {params.extra} 1> {log} 2>&1
        """

rule make_simulation_transcriptome:
    input:
        fasta=config["genome"],
        annotation="results/filter_annotation/filtered_annotation.gtf",
    output:
        records="results/make_simulation_transcriptome/transcriptome_sim.fasta",
    threads: 1
    log:
        "logs/make_simulation_transcriptome/gffread.log",
    params:
        extra=config["simulation_gffread_extra"]
    conda:
        "../envs/gffread.yaml"
    script: 
        "../scripts/gffread.py"



rule simulate_fastas:
    input:
        fasta="results/make_simulation_transcriptome/transcriptome_sim.fasta",
        counts="results/filter_annotation/transcript_read_counts.csv",
    output:
        sim=SIMULATION_OUTPUT
    threads: 1
    log:
        "logs/simulate_fastas/simulate_fastas.log"
    conda:
        "../envs/simulate.yaml"
    params:
        rscript=workflow.source_path("../scripts/simulate.R"),
        nreps=config["number_of_replicates"],
        extra=SIMULATION_PARAMS
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -f {input.fasta} -c {input.counts} -o ./results/simulate_fastas/ \
        -n {params.nreps} {params.extra} 1> {log} 2>&1
        """

