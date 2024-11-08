"""
Annotate TnpB protein using pfam domains and identify conserved IS region.
"""

import textwrap


rule hmmer_search:
    input:
        fasta="results/genomes/{genome}/translated_cds.faa",
        profile=config["domain_search"]["models"],
    output:
        domtblout="results/genomes/{genome}/hmmer.txt",
    log:
        hmmer="results/genomes/{genome}/hmmer.log",
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
        hits="results/genomes/{genome}/hmmer_hits.bed",
    conda:
        "../env/csvtk_gawk.yaml"
    shell:
        textwrap.dedent(r"""
        cat {input.hmmer:q} |
            gawk '!/^#/ && $7<0.001 && $12<0.001 {{
                match($0, /\[locus_tag=([^[:space:]]+)\]/, arr)
                print arr[1] "\t" $4
            }}' |
            {{ csvtk -tH fold -f 1 -v 2 || echo '#' ; }} |
            csvtk -tH join --left-join - {input.annotation:q} -f '1;12' |
            awk -v 'FS=\t' -v 'OFS=\t' '{{ print $4, $6, $7, $1, ".", $5, $2 }}' |
            sort -k1,1 -k2,2n > {output.hits:q}
        """)


rule match_tnp:
    input:
        genes=rules.hmmer_filter.output,
        fasta="results/genomes/{genome}/translated_cds.faa",
    output:
        bed="results/genomes/{genome}/{tnp}.bed",
        fasta="results/genomes/{genome}/{tnp}.faa",
    wildcard_constraints:
        tnp="tnp[ab]"
    params:
        motifs=lambda wildcards: config["domain_search"]["rule"][wildcards.tnp],
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

