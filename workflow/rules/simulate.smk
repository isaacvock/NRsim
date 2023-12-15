rule filter_annotation:
    input:
        gtf="results/modify_annotation/modified_annotation.gtf"
    output:
        filtered="results/filter_annotation/filtered_annotation.gtf"
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
        {params.rscript} -g {input.gtf} -o {output.filtered} {params.extra} 1> {log} 2>&1
        """




rule simulate_fastas:



rule convert_to_fastq: