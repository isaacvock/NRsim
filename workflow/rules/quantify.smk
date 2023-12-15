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
    log:
        "logs/modify_annotation/modify_annotation.log",
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


### Quantify with Salmon

# Never make decoy index because I want to quantify intronic content accurately
rule index:
    input:
        sequences="results/make_transcriptome_fasta/transcriptome.fasta",
    output:
        multiext(
            INDEX_PATH,
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
            r1="results/trimmed/read.1.fastq",
            r2="results/trimmed/read.2.fastq",
            index=multiext(
                    INDEX_PATH,
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
            quant="results/quant/quant.sf",
            lib="results/quant/lib_format_counts.json"
        log:
            "logs/quant/salmon.log"
        params:
            libtype=LIBTYPE,
            extra=config["salmon_quant_params"]
        threads: 12 # See https://salmon.readthedocs.io/en/latest/salmon.html Note for motivation
        wrapper:
            "v2.6.0/bio/salmon/quant"

else:

    rule quant:
        input:
            r="results/trimmed/read.1.fastq",
            index=multiext(
                    INDEX_PATH,
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
            quant="results/quant/quant.sf",
            lib="results/quant/lib_format_counts.json"
        log:
            "logs/quant/salmon.log"
        params:
            libtype=LIBTYPE,
            extra=config["salmon_quant_params"]
        threads: 12
        wrapper:
            "v2.6.0/bio/salmon/quant"
