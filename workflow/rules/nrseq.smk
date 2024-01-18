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

        rule split_fastas:
            input:
                fasta="results/simulate_fastas/{sim}/sample_{sample}_{read}.fasta",
            output:
                temp(expand("results/split_fastas/{{sim}}/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/{sim}/sample_{sample}_read_{read}.log"
            params:
                nsplit = NUMBER_SPLIT,
                reads = lambda w: config["simulation_parameters"]["library_size"]["{}".format(w.sim)],
                outdir = lambda w: "results/split_fastas/{}/".format(w.sim),
            conda:
                "../envs/fastq.yaml"
            threads: 1
            script:
                "../scripts/split_fasta.py"

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/split_fastas/{{sim}}/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS),
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

        rule shuffle_fastq:
            input:
                expand("results/convert_to_fastq/{{sim}}/sample_{{sample}}_{READ}.fastq.gz", READ = ['1', '2']),
            output:
                tread1 = temp("results/shuffle_fastq/{sim}/sample_{sample}/read_1.fastq"),
                tread2 = temp("results/shuffle_fastq/{sim}/sample_{sample}/read_2.fastq"),
                read1 = "results/shuffle_fastq/{sim}/sample_{sample}/read_1.fastq.gz",
                read2 = "results/shuffle_fastq/{sim}/sample_{sample}/read_2.fastq.gz",
            log:
                "logs/shuffle_fastq/{sim}/sample_{sample}.log",
            conda:
                "../envs/fastq.yaml"
            params:
                shellscript = workflow.source_path("../scripts/shuffle_pe.sh"),
            threads: 8
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {input} {output.tread1} {output.tread2} {output.read1} {output.read2} 1> {log} 2>&1
                """

    else:

        rule split_fastas:
            input:
                fasta="results/simulate_fastas/{sim}/sample_{sample}.fasta",
            output:
                temp(expand("results/split_fastas/{{sim}}/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/{sim}/sample_{sample}.log"
            params:
                nsplit = NUMBER_SPLIT,
                reads = lambda w: config["simulation_parameters"]["library_size"]["{}".format(w.sim)],
                outdir = lambda w: "results/split_fastas/{}/".format(w.sim),
            conda:
                "../envs/fastq.yaml"
            threads: 1
            script:
                "../scripts/split_fasta.py"

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/split_fastas/{{sim}}/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS),
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
                {output.fastq} {params.PE} ./results/simulate_fastas/{wildcards.sim}/ {params.mutrate} DUMMY 1> {log} 2>&1
                """

        rule shuffle_fastq:
            input:
                "results/convert_to_fastq/{sim}/sample_{sample}.fastq.gz",
            output:
                tread = temp("results/shuffle_fastq/{sim}/sample_{sample}/read_1.fastq"),
                read = "results/shuffle_fastq/{sim}/sample_{sample}/read_1.fastq.gz",
            log:
                "logs/shuffle_fastq/{sim}/sample_{sample}.log"
            conda:
                "../envs/fastq.yaml"
            params:
                shellscript = workflow.source_path("../scripts/shuffle_se.sh"),
            threads: 8
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {input} {output.tread} {output.read} {threads} 1> {log} 2>&1
                """

else:

    if PE:

        rule split_fastas:
            input:
                fasta="results/simulate_fastas/sample_{sample}_{read}.fasta",
            output:
                temp(expand("results/split_fastas/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/sample_{sample}_read_{read}.log"
            params:
                nsplit = NUMBER_SPLIT,
                reads = config["library_size"],
                outdir = lambda w: "results/split_fastas/{}/".format(w.sim),
            conda:
                "../envs/fastq.yaml"
            threads: 1
            script:
                "../scripts/split_fasta.py"

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/split_fastas/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS),
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

        rule shuffle_fastq:
            input:
                expand("results/convert_to_fastq/sample_{{sample}}_{READ}.fastq.gz", READ = ['1', '2']),
            output:
                tread1 = temp("results/shuffle_fastq/sample_{sample}/read_1.fastq"),
                tread2 = temp("results/shuffle_fastq/sample_{sample}/read_2.fastq"),
                read1 = "results/shuffle_fastq/sample_{sample}/read_1.fastq.gz",
                read2 = "results/shuffle_fastq/sample_{sample}/read_2.fastq.gz",
            log:
                "logs/shuffle_fastq/sample_{sample}.log",
            conda:
                "../envs/fastq.yaml"
            params:
                shellscript = workflow.source_path("../scripts/shuffle_pe.sh"),
            threads: 8
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {input} {output.tread1} {output.tread2} {output.read1} {output.read2} 1> {log} 2>&1
                """


    else:

        rule split_fastas:
            input:
                fasta="results/simulate_fastas/sample_{sample}.fasta",
            output:
                temp(expand("results/split_fastas/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS)),
            log:
                "logs/split_fasta/sample_{sample}.log"
            params:
                nsplit = NUMBER_SPLIT,
                reads = config["library_size"],
                outdir = lambda w: "results/split_fastas/{}/".format(w.sim),
            conda:
                "../envs/fastq.yaml"
            threads: 1
            script:
                "../scripts/split_fasta.py"

        rule convert_to_fastq:
            input:
                kinetics="results/generate_transcript_kinetics/kinetics.csv",
                fasta=expand("results/split_fastas/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS),
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

        rule shuffle_fastq:
            input:
                "results/convert_to_fastq/sample_{sample}.fastq.gz",
            output:
                tread = temp("results/shuffle_fastq/sample_{sample}/read_1.fastq"),
                read = "results/shuffle_fastq/sample_{sample}/read_1.fastq.gz",
            log:
                "logs/shuffle_fastq/sample_{sample}.log"
            conda:
                "../envs/fastq.yaml"
            params:
                shellscript = workflow.source_path("../scripts/shuffle_se.sh"),
            threads: 8
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {input} {output.tread} {output.read} {threads} 1> {log} 2>&1
                """
