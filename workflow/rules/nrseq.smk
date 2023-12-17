if PE:

    rule split_fasta:
        input:
            fasta="results/simulate_fastas/sample_{{sample}}_{{read}}.fasta",
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
            pyfasta split -n {params.nsub} {input.fasta}
            """

    rule convert_to_fastq:
        input:
            fasta=expand("results/simulate_fastas/sample_{{sample}}_{{read}}.{ID}.fasta", ID = SPLIT_IDS),
        output:
            fastq="results/convert_to_fastq/sample_{sample}_{read}.fastq.gz"
        log:
            "logs/convert_to_fastq/sample_{sample}_{read}.log"
        params:
            qscore=config["qscore_impute"],
            shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
            pythonscript=workflow.source_path("../scripts/fasta_to_fastq.py")
        conda:
            "../envs/fastq.yaml"
        threads: 32
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} {output.fastq} {wildcards.read} 1> {log} 2>&1
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
            pyfasta split -n {params.nsub} {input.fasta}
            """

    rule convert_to_fastq:
        input:
            fasta=expand("results/simulate_fastas/sample_{{sample}}.{ID}.fasta", ID = SPLIT_IDS),
        output:
            fastq="results/convert_to_fastq/sample_{sample}.fastq.gz"
        log:
            "logs/convert_to_fastq/sample_{sample}.log"
        params:
            qscore=config["qscore_impute"],
            shellscript=workflow.source_path("../scripts/fasta_to_fastq.sh"),
            pythonscript=workflow.source_path("../scripts/fasta_to_fastq.py")
        conda:
            "../envs/fastq.yaml"
        threads: 32
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {params.qscore} {params.pythonscript} {output.fastq} DUMMY 1> {log} 2>&1
            """


#rule generate_transcript_kinetics:


#rule make_nrseq_fastq: