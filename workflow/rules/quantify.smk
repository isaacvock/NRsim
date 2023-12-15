# Make modified gtf for qunatification
rule modify_annotation:
    input:
        gtf=config["annotation"]
    output:
        mod_gtf="results/modify_annotation/modified_annotation.gtf"
    params:
        rscript=workflow.source_path("../scripts/modify_annotation.R"),
        extra=config["modify_annotation_params"]
    threads: 1
    conda:
        "../envs/simulate.yaml"
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -g {input.gtf} -o {output.mod_gtf} {params.extra}
        """

# Make transcriptome fasta from genome fasta + gtf
rule make_transcriptome_fasta:
    input:
        fasta=config["genome"],
        annotation="results/modify_gtf/modified_annotation.gtf",
    output:
        records="results/make_transcriptome_fasta/transcriptome.fasta",
    threads: 1
    log:
        "logs/make_transcriptome_fasta/gffread.log",
    params:
        extra=config["gffread_extra"]
    conda:
        "../envs/gffread.yaml"
    script: 
        "../scripts/gffread.py"


if config["quantifier"] == "kallisto":

    ### Quantify with Kallisto


    rule index:
        input:
            fasta="results/make_transcriptome_fasta/transcriptome.fasta",
        output:
            index=INDEX_KALLISTO,
        params:
            extra=config["kallisto_index_params"]
        log:
            "logs/index/kallisto_index.log"
        threads: 8
        wrapper:
            "v2.6.0/bio/kallisto/index"

    rule quant:
        input:
            fastq=expand("results/trimmed/{{sample}}.{read}.fastq", read = READS),
            index=INDEX_KALLISTO
        output:
            dir=directory("results/kallisto_quant/{sample}"),
            run_info="results/kallisto_quant/{sample}/run_info.json"
        params:
            extra=config["kallisto_quant_params"],
        log:
            "logs/quant/{sample}_kallisto.log"
        conda:
            "../envs/kallisto.yaml"
        threads: 12
        script:
            "../scripts/kallisto-quant.py"
        

else:

    ### Quantify with Salmon

    # Never make decoy index because I want to quantify intronic content accurately
    rule index:
        input:
            sequences=SALMON_TRANSCRIPTOME,
        output:
            multiext(
                config["indices"],
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        log:
            "logs/index/salmon_index.log"
        threads: 16 # Borrowed from vignette here: https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
        params:
            extra=config["salmon_index_params"]
        wrapper:
            "v2.6.0/bio/salmon/index"



    if config["PE"]:

        rule quant:
            input:
                r1=get_fastq_r1,
                r2=get_fastq_r2,
                index=multiext(
                        config["indices"],
                        "complete_ref_lens.bin",
                        "ctable.bin",
                        "ctg_offsets.bin",
                        "duplicate_clusters.tsv",
                        "info.json",
                        "mphf.bin",
                        "pos.bin",
                        "pre_indexing.log",
                        "rank.bin",
                        "refAccumLengths.bin",
                        "ref_indexing.log",
                        "reflengths.bin",
                        "refseq.bin",
                        "seq.bin",
                        "versionInfo.json",
                    ),
            output:
                quant="results/quant/{sample}/quant.sf",
                lib="results/quant/{sample}/lib_format_counts.json"
            log:
                "logs/quant/{sample}_salmon.log"
            params:
                libtype=LIBTYPE,
                extra=config["salmon_quant_params"]
            threads: 12 # See https://salmon.readthedocs.io/en/latest/salmon.html Note for motivation
            wrapper:
                "v2.6.0/bio/salmon/quant"

    else:

        rule quant:
            input:
                r=get_fastq_r1,
                index=multiext(
                        config["indices"],
                        "complete_ref_lens.bin",
                        "ctable.bin",
                        "ctg_offsets.bin",
                        "duplicate_clusters.tsv",
                        "info.json",
                        "mphf.bin",
                        "pos.bin",
                        "pre_indexing.log",
                        "rank.bin",
                        "refAccumLengths.bin",
                        "ref_indexing.log",
                        "reflengths.bin",
                        "refseq.bin",
                        "seq.bin",
                        "versionInfo.json",
                    ),
            output:
                quant="results/quant/{sample}/quant.sf",
                lib="results/quant/{sample}/lib_format_counts.json"
            log:
                "logs/quant/{sample}_salmon.log"
            params:
                libtype=LIBTYPE,
                extra=config["salmon_quant_params"]
            threads: 12
            wrapper:
                "v2.6.0/bio/salmon/quant"
