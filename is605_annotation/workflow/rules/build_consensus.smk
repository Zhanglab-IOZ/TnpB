"""
Build consensus from multiple IS copies.
"""

import textwrap


rule blast:
    input:
        copies="results/elements_list/{tnpb}/flanking/{tnpb}.fna",
    output:
        table="results/elements/{tnpb}/blast.tsv",
    conda: "../env/blast.yaml"
    shell:
        "blastn -query {input.copies:q} -subject {input.copies:q} -outfmt '7 qaccver qlen saccver slen pident length qstart qend sstart send evalue bitscore ' > {output:q}"


rule clip:
    input:
        blast=rules.blast.output.table,
        copies="results/elements_list/{tnpb}/flanking/{tnpb}.fna",
    output:
        bed="results/elements/{tnpb}/is.bed",
        fasta="results/elements/{tnpb}/is.fna",
    conda: "../env/bedtools_csvtk.yaml"
    shell:
        textwrap.dedent(r"""
            cat {input.blast:q} |
                csvtk -tH filter2 -f ' $7 - 1 < 1500 && $2 - $8 < 1500 && $6 < 2500 ' |
                {{ csvtk -tH summary -g 1 -f '7:min,8:max' -n 0 || true; }} |
                tee {output.bed:q} |
                bedtools getfasta -s -fi {input.copies:q} -bed - -bedOut | 
                awk -v 'FS=\t' '{{ print ">"$1 "\n" $4 }}' >{output.fasta:q}
        """)


rule mafft:
    input:
        fasta=rules.clip.output.fasta,
    output:
        fasta="results/elements/{tnpb}/mafft.fna",
    conda: "../env/mafft.yaml"
    log:
        mafft = "results/elements/{tnpb}/mafft.log",
    shell:
        "mafft-linsi {input.fasta:q} > {output.fasta:q} 2>{log.mafft:q}"
