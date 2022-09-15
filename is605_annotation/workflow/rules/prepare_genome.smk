import textwrap


rule prepare_genome:
    input:
        archive=lambda wildcards: config["input"]["genomes"][wildcards.genome].format(suffix = wildcards.suffix),
    output:
        extracted=temp("results/genomes/{genome}/{suffix}"),
    wildcard_constraints:
        suffix="translated_cds.faa|genomic.gff|genomic.fna"
    shell:
        " zcat {input.archive:q} > {output.extracted:q} "


rule prepare_annotation:
    input:
        gff = "results/genomes/{genome}/genomic.gff",
    output:
        genepred = temp("results/genomes/{genome}/annotation.genepred"),
    conda:
        "../env/gff3togenepred.yaml"
    shell:
        "gff3ToGenePred {input.gff:q} {output.genepred:q} -geneNameAttr=locus_tag -useName"


rule prepare_genome_index:
    input:
        genome="results/genomes/{genome}/genomic.fna",
    output:
        index= temp("results/genomes/{genome}/genomic.fna.fai"),
    conda:
        "../env/seqkit.yaml"
    shell:
        "seqkit faidx {input.genome:q}"


