rule filter_annotation:
    input:
        gtf="results/modify_annotation/modified_annotation.gtf",
        quantification="results/quant/quant.sf",
    output:
        filtered="results/filter_annotation/filtered_annotation.gtf",
        counts="results/filter_annotation/transcript_read_counts.csv",
        nointron="results/filter_annotation/filtered_annotation_nointron.gtf",
    params:
        rscript=workflow.source_path("../scripts/filter_annotation.R"),
        extra=config["filter_annotation_params"],
    conda:
        "../envs/simulate.yaml"
    log:
        "logs/filter_annotation/filter_annotation.log",
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -g {input.gtf} -q {input.quantification} \
        -c {output.counts} -o {output.filtered} -n {output.nointron} {params.extra} 1> {log} 2>&1
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
        extra=config["simulation_gffread_extra"],
    conda:
        "../envs/gffread.yaml"
    script:
        "../scripts/gffread.py"


if config["dataset_specific"]:
    if PE:

        rule simulate_fastas:
            input:
                fasta="results/make_simulation_transcriptome/transcriptome_sim.fasta",
                counts="results/filter_annotation/transcript_read_counts.csv",
            output:
                expand(
                    "results/simulate_fastas/{{sim}}/sample_{SID}_{READS}.fasta",
                    SID=sample_names,
                    READS=READS,
                ),
            threads: 1
            log:
                "logs/simulate_fastas/{sim}/simulate_fastas.log",
            conda:
                "../envs/simulate.yaml"
            params:
                rscript=workflow.source_path("../scripts/simulate.R"),
                nreps=config["number_of_replicates"],
                library_size=lambda wildcards: config["simulation_parameters"][
                    "library_size"
                ][str(wildcards.sim)],
                sim_premrna=lambda wildcards: config["simulation_parameters"][
                    "simulate_premRNA"
                ][str(wildcards.sim)],
                extra=lambda wildcards: config["simulation_parameters"][
                    "extra_params"
                ][str(wildcards.sim)],
            shell:
                """
                chmod +x {params.rscript}
                {params.rscript} -f {input.fasta} -c {input.counts} -o ./results/simulate_fastas/{wildcards.sim} \
                -n {params.nreps} -l {params.library_size} -m {params.sim_premrna} {params.extra} 1> {log} 2>&1
                """

    else:

        rule simulate_fastas:
            input:
                fasta="results/make_simulation_transcriptome/transcriptome_sim.fasta",
                counts="results/filter_annotation/transcript_read_counts.csv",
            output:
                expand(
                    "results/simulate_fastas/{{sim}}/sample_{SID}.fasta",
                    SID=sample_names,
                ),
            threads: 1
            log:
                "logs/simulate_fastas/{sim}/simulate_fastas.log",
            conda:
                "../envs/simulate.yaml"
            params:
                rscript=workflow.source_path("../scripts/simulate.R"),
                nreps=config["number_of_replicates"],
                library_size=lambda wildcards: config["simulation_parameters"][
                    "library_size"
                ][str(wildcards.sim)],
                sim_premrna=lambda wildcards: config["simulation_parameters"][
                    "simulate_premRNA"
                ][str(wildcards.sim)],
                extra=lambda wildcards: config["simulation_parameters"][
                    "extra_params"
                ][str(wildcards.sim)],
            shell:
                """
                chmod +x {params.rscript}
                {params.rscript} -f {input.fasta} -c {input.counts} -o ./results/simulate_fastas/{wildcards.sim} \
                -n {params.nreps} -l {params.library_size} -m {params.sim_premrna} --singleend {params.extra} 1> {log} 2>&1
                """

else:

    rule simulate_fastas:
        input:
            fasta="results/make_simulation_transcriptome/transcriptome_sim.fasta",
            counts="results/filter_annotation/transcript_read_counts.csv",
        output:
            sim=SIMULATION_OUTPUT,
        threads: 1
        log:
            "logs/simulate_fastas/simulate_fastas.log",
        conda:
            "../envs/simulate.yaml"
        params:
            rscript=workflow.source_path("../scripts/simulate.R"),
            nreps=config["number_of_replicates"],
            library_size=config["library_size"],
            sim_premrna=config["simulate_premRNA"],
            pe=lambda wildcards: "" if config["PE"] else "--singleend",
            extra=config["simulate_fastas_params"],
        shell:
            """
            chmod +x {params.rscript}
            {params.rscript} -f {input.fasta} -c {input.counts} -o ./results/simulate_fastas/ \
            -n {params.nreps} -l {params.library_size} -m {params.sim_premrna} {params.pe} {params.extra} 1> {log} 2>&1
            """
