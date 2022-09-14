"""
Annotate TnpB protein using pfam domains and identify conserved IS region.
"""

import textwrap


rule prepare_genome:
    input:
        archive=lambda wildcards: config["input"]["genomes"][wildcards.genome].format(suffix = wildcards.suffix),
    output:
        extracted=temp("results/{genome}/{suffix}"),
    wildcard_constraints:
        suffix="translated_cds.faa|genomic.gff|genomic.fna"
    shell:
        " zcat {input.archive:q} > {output.extracted:q} "


rule prepare_annotation:
    input:
        gff = "results/{genome}/genomic.gff",
    output:
        genepred = temp("results/{genome}/annotation.genepred"),
    conda:
        "../env/gff3togenepred.yaml"
    shell:
        "gff3ToGenePred {input.gff:q} {output.genepred:q} -geneNameAttr=locus_tag -useName"


rule prepare_genome_index:
    input:
        genome="results/{genome}/genomic.fna",
    output:
        index= temp("results/{genome}/genomic.fna.fai"),
    conda:
        "../env/seqkit.yaml"
    shell:
        "seqkit faidx {input.genome:q}"


rule hmmer_search:
    input:
        fasta="results/{genome}/translated_cds.faa",
        profile=config["params"]["pfam"]["models"],
    output:
        domtblout="results/{genome}/hmmer.txt",
    log:
        hmmer="results/{genome}/hmmer.log",
    conda:
        "../env/hmmer.yaml"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --domtblout {output.domtblout:q} {input.profile:q} {input.fasta:q} | gzip > {log.hmmer:q}"


rule hmmer_filter:
    input:
        hmmer=rules.hmmer_search.output.domtblout,
        annotation=rules.prepare_annotation.output,
    output:
        hits="results/{genome}/hmmer_hits.bed",
    conda:
        "../env/csvtk.yaml"
    shell:
        textwrap.dedent(r"""
        cat {input.hmmer:q} |
            awk '!/^#/ && $7<0.001 && $12<0.001 {{
                match($0, /\[locus_tag=([^[:space:]]+)\]/, arr)
                print arr[1] "\t" $4
            }}' |
            csvtk -tH fold -f 1 -v 2 |
            csvtk -tH join --left-join - {input.annotation:q} -f '1;12' |
            awk -v 'FS=\t' -v 'OFS=\t' '{{ print $4, $6, $7, $1, ".", $5, $2 }}' |
            sort -k1,1 -k2,2n > {output.hits:q}
        """)


rule match_tnp:
    input:
        genes=rules.hmmer_filter.output,
        fasta="results/{genome}/translated_cds.faa",
    output:
        bed="results/{genome}/{tnp}.bed",
        fasta="results/{genome}/{tnp}.faa",
    params:
        motifs=lambda wildcards: config["params"]["pfam"]["rule"][wildcards.tnp],
    conda:
        "../env/seqkit.yaml"
    shell:
        textwrap.dedent(r"""
            awk -v 'FS=\t' '$7 ~ /{params.motifs}/' {input.genes:q} |
                tee {output.bed:q} |
                cut -f 4 |
                seqkit grep -f - --id-regexp '\[locus_tag=(\S+)\]' {input.fasta:q} |
                seqkit seq -i --id-regexp '\[locus_tag=(\S+)\]' > {output.fasta:q}
        """)

