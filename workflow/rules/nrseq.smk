rule generate_transcript_kinetics:
    input:
        counts="results/filter_annotation/transcript_read_counts.csv",
    output:
        kinetics="results/generate_transcript_kinetics/kinetics.csv"
    log:
        "logs/generate_transcript_kinetics/generate_transcript_kinetics.log"
    conda:
        "../envs/simulate.yaml"
    params:
        rscript=workflow.source_path("../scripts/generate_transcript_kinetics.R"),
        extra=config["generate_transcript_kinetics_params"]
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -c {input.counts} -o {output.kinetics} {params.extra} 1> {log} 2>&1
        """

if config["dataset_specific"]:

    if PE:

        rule split_fasta:
            input:
                fasta="results/simulate_fastas/{sim}/sample_{sample}_{read}.fasta",
            output:
                temp(expand("results/simulate_fastas/{{sim}}/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/{sim}/sample_{sample}_read_{read}.log"
            params:
                nsub=NUMBER_SPLIT,
            conda:
                "../envs/pyfasta.yaml"
            threads: 1
            shell:
                """
                pyfasta split -n {params.nsub} {input.fasta}
                """

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/simulate_fastas/{{sim}}/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS),
            output:
                fastq="results/convert_to_fastq/{sim}/sample_{sample}_{read}.fastq.gz",
            log:
                "logs/convert_to_fastq/{sim}/sample_{sample}_{read}.log"
            params:
                qscore=lambda wildcards: config["simulation_parameters"]["qscore_impute"][str(wildcards.sim)],
                shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
                pythonscript=workflow.source_path("../scripts/make_nrseq_fastq.py"),
                PE = PE,
                mutrate = lambda wildcards: config["simulation_parameters"]["s4U_mutation_rate"][str(wildcards.sim)]
            conda:
                "../envs/fastq.yaml"
            threads: 32
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} \
                {input.kinetics} {output.fastq} {params.PE} ./results/simulate_fastas/{wildcards.sim}/ {params.mutrate} {wildcards.read} 1> {log} 2>&1
                """

    else:

        rule split_fasta:
            input:
                fasta="results/simulate_fastas/{sim}/sample_{sample}.fasta",
            output:
                expand("results/simulate_fastas/{{sim}}/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS),
            log:
                "logs/split_fasta/{sim}/sample_{sample}.log"
            params:
                nsub=NUMBER_SPLIT,
            conda:
                "../envs/pyfasta.yaml"
            threads: 1
            shell:
                """
                pyfasta split -n {params.nsub} {input.fasta}  1> {log} 2>&1
                """

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=temp(expand("results/simulate_fastas/{{sim}}/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS)),
            output:
                fastq="results/convert_to_fastq/{sim}/sample_{sample}.fastq.gz",
            log:
                "logs/convert_to_fastq/{sim}/sample_{sample}.log"
            params:
                qscore= lambda wildcards: config["simulation_parameters"]["qscore_impute"][str(wildcards.sim)],
                shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
                pythonscript=workflow.source_path("../scripts/make_nrseq_fastq.py"),
                PE = PE,
                mutrate = lambda wildcards: config["simulation_parameters"]["s4U_mutation_rate"][str(wildcards.sim)]
            conda:
                "../envs/fastq.yaml"
            threads: 32
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} {input.kinetics} \
                {output.fastq} {params.PE} ./results/simulate_fastas/{wildcards.sim}/ {params.mutrate} {wildcards.read} 1> {log} 2>&1
                """

else:

    if PE:

        rule split_fasta:
            input:
                fasta="results/simulate_fastas/sample_{sample}_{read}.fasta",
            output:
                temp(expand("results/simulate_fastas/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/sample_{sample}_read_{read}.log"
            params:
                nsub=NUMBER_SPLIT,
            conda:
                "../envs/pyfasta.yaml"
            threads: 1
            shell:
                """
                pyfasta split -n {params.nsub} {input.fasta}  1> {log} 2>&1
                """

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/simulate_fastas/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS),
            output:
                fastq="results/convert_to_fastq/sample_{sample}_{read}.fastq.gz",
            log:
                "logs/convert_to_fastq/sample_{sample}_{read}.log"
            params:
                qscore=config["qscore_impute"],
                shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
                pythonscript=workflow.source_path("../scripts/make_nrseq_fastq.py"),
                PE = PE,
                mutrate = config["s4U_mutation_rate"]
            conda:
                "../envs/fastq.yaml"
            threads: 32
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} {input.kinetics} \
                {output.fastq} {params.PE} ./results/simulate_fastas/ {params.mutrate} {wildcards.read} 1> {log} 2>&1
                """


    else:

        rule split_fasta:
            input:
                fasta="results/simulate_fastas/sample_{sample}.fasta",
            output:
                temp(expand("results/simulate_fastas/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/sample_{sample}.log"
            params:
                nsub=NUMBER_SPLIT,
            conda:
                "../envs/pyfasta.yaml"
            threads: 1
            shell:
                """
                pyfasta split -n {params.nsub} {input.fasta}  1> {log} 2>&1
                """

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/simulate_fastas/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS),
            output:
                fastq="results/convert_to_fastq/sample_{sample}.fastq.gz"
            log:
                "logs/convert_to_fastq/sample_{sample}.log"
            params:
                qscore=config["qscore_impute"],
                shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
                pythonscript=workflow.source_path("../scripts/make_nrseq_fastq.py"),
                PE = PE,
                mutrate = config["s4U_mutation_rate"]
            conda:
                "../envs/fastq.yaml"
            threads: 32
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} \
                {input.kinetics} {output.fastq} {params.PE} ./results/simulate_fastas/ {params.mutrate} DUMMY 1> {log} 2>&1
                """