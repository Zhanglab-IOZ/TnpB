"""
Build consensus from multiple IS copies.
"""

import textwrap


rule blast:
    input:
        copies="results/flanking/{tnpb}/flanking.fna",
    output:
        table="results/flanking/{tnpb}/blast.tsv",
    conda: "../env/blast.yaml"
    shell:
        "blastn -query {input.copies:q} -subject {input.copies:q} -outfmt '7 qaccver qlen saccver slen pident length qstart qend sstart send evalue bitscore ' > {output:q}"


rule clip:
    input:
        blast=rules.blast.output.table,
        copies="results/flanking/{tnpb}/flanking.fna",
    output:
        bed="results/flanking/{tnpb}/is.bed",
        fasta="results/flanking/{tnpb}/is.fna",
    conda: "../env/bedtools_csvtk.yaml"
    shell:
        textwrap.dedent(r"""
            cat {input.blast:q} |
                csvtk -tH filter2 -f ' $7 - 1 < 1500 && $2 - $8 < 1500 && $6 < 2500 ' |
                csvtk -tH summary -g 1 -f '7:min,8:max' -n 0 |
                tee {output.bed:q} |
                bedtools getfasta -s -fi {input.copies:q} -bed - -bedOut | 
                awk -v 'FS=\t' '{{ print ">"$1 "\n" $4 }}' >{output.fasta:q}
        """)


rule mafft:
    input:
        fasta=rules.clip.output.fasta,
    output:
        fasta="results/flanking/{tnpb}/mafft.fna",
    conda: "../env/mafft.yaml"
    log:
        mafft = "results/flanking/{tnpb}/mafft.log",
    shell:
        "mafft-linsi {input.fasta:q} > {output.fasta:q} 2>{log.mafft:q}"
