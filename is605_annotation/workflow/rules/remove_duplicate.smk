import textwrap


rule merge_tnpb:
    input:
        sequences = expand("results/genomes/{genome}/multi_copy.faa", genome = config["input"]["genomes"].keys()),
    output:
        "results/dedup/denovo_tnpb.faa",
    shell:
        r"""
            cat {input:q} >{output:q}
        """

