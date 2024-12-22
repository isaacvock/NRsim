# Unzip if gzipped
if is_gz:

    rule unzip:
        input:
            fastqs=FASTQ_FILES,
        output:
            unzipped_fqs=temp(
                expand("results/unzipped/read_{read}.fastq", read=READS_input)
            ),
        log:
            "logs/unzip/unzip.log",
        conda:
            "../envs/pigz.yaml"
        threads: 10
        script:
            "../scripts/pigz.py"


## Trim adapters
if PE_input:
    if is_gz:

        # Trim with fastp (automatically detects adapters)
        rule fastp:
            input:
                sample=expand("results/unzipped/read_{read}.fastq", read=READS_input),
            output:
                trimmed=["results/trimmed/read.1.fastq", "results/trimmed/read.2.fastq"],
                # Unpaired reads separately
                unpaired1="results/trimmed/read.u1.fastq",
                unpaired2="results/trimmed/read.u2.fastq",
                failed="results/trimmed/reads.failed.fastq",
                html="results/reports/reads.html",
                json="results/reports/reads.json",
            log:
                "logs/fastp/fastp.log",
            params:
                adapters=config["fastp_adapters"],
                extra="",
            threads: 20
            wrapper:
                "v2.2.1/bio/fastp"

    else:

        # Trim with fastp (automatically detects adapters)
        rule fastp:
            input:
                sample=FASTQ_FILES,
            output:
                trimmed=["results/trimmed/read.1.fastq", "results/trimmed/read.2.fastq"],
                # Unpaired reads separately
                unpaired1="results/trimmed/read.u1.fastq",
                unpaired2="results/trimmed/read.u2.fastq",
                failed="results/trimmed/reads.failed.fastq",
                html="results/reports/reads.html",
                json="results/reports/reads.json",
            log:
                "logs/fastp/fastp.log",
            params:
                adapters=config["fastp_adapters"],
                extra="",
            threads: 20
            wrapper:
                "v2.2.1/bio/fastp"

else:
    if is_gz:

        # Trim with fastp (automatically detects adapters)
        rule fastp:
            input:
                sample=expand("results/unzipped/read_{read}.fastq", read=READS_input),
            output:
                trimmed="results/trimmed/read.1.fastq",
                failed="results/trimmed/read.1.failed.fastq",
                html="results/reports/read.1.html",
                json="results/reports/read.1.json",
            log:
                "logs/fastp/fastp.log",
            params:
                adapters=config["fastp_adapters"],
                extra="",
            threads: 20
            wrapper:
                "v2.2.1/bio/fastp"

    else:

        # Trim with fastp (automatically detects adapters)
        rule fastp:
            input:
                sample=FASTQ_FILES,
            output:
                trimmed="results/trimmed/read.1.fastq",
                failed="results/trimmed/read.1.failed.fastq",
                html="results/reports/read.1.html",
                json="results/reports/read.1.json",
            log:
                "logs/fastp/fastp.log",
            params:
                adapters=config["fastp_adapters"],
                extra="",
            threads: 20
            wrapper:
                "v2.2.1/bio/fastp"


# Run fastqc on trimmed fastqs
rule fastqc:
    input:
        "results/trimmed/read.{read}.fastq",
    output:
        html="results/fastqc/read{read}.html",
        zip="results/fastqc/read{read}_fastqc.zip",
    log:
        "logs/fastqc/read{read}.log",
    params:
        extra=config["fastqc_params"],
    resources:
        mem_mb=9000,
    threads: 4
    wrapper:
        "v2.2.1/bio/fastqc"
