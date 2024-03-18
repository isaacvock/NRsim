

rule wtf:
    output:
        "wtf_{sim}_{sample}_{read}.txt"
    params:
        seed=lambda wildcards: seeds[str(wildcards.sim)][str(wildcards.sample)],
    shell:
        """
        echo {params.seed} > {output}
        """