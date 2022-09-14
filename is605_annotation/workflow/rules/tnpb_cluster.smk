import textwrap


rule cluster_tnpb:
    input:
        fasta="results/{genome}/tnpb.faa",
    output:
        rep="results/{genome}/tnpb_rep_seq.fasta",
        allseqs="results/{genome}/tnpb_all_seqs.fasta",
        table="results/{genome}/tnpb_cluster.tsv",
    log:
        mmseqs2="results/{genome}/tnpb_mmseqs.log",
    params:
        prefix="results/{genome}/tnpb",
        identity=0.9,
        coverage=0.9,
    conda:
        "../env/mmseqs2.yaml"
    threads: 8
    shell:
        textwrap.dedent(r"""
            tmpd=$(mktemp -d -t --dry-run mmseq_XXXXXXXX)
            mmseqs easy-cluster {input:q} {params.prefix:q} $tmpd --min-seq-id {params.identity} -c {params.coverage} --threads {threads} > {log.mmseqs2:q}
            rm -r $tmpd
        """)


rule summary_tnpb:
    input:
        tnpa="results/{genome}/tnpa.bed",
        tnpb="results/{genome}/tnpb.bed",
        tnpb_cluster=rules.cluster_tnpb.output.table,
    output:
        details="results/{genome}/tnpb_annot.bed",
        summary="results/{genome}/summary.tsv",
    conda:
        "../env/bedtools_csvtk.yaml"
    shell:
        textwrap.dedent(r"""
            bedtools closest -a {input.tnpb:q} -b {input.tnpa:q} -d |
                cut -f 1-7,15 |
                awk -v 'FS=\t' -v 'OFS=\t' '{{ $8 = NF>=8 ? ($8<500 && $8!=-1) : 0 ; print }}' |
                csvtk -tH join --left-join - {input.tnpb_cluster:q} -f "4;2" |
                tee {output.details:q} |
                csvtk -tH summary -g 9 -f '4:count,4:collapse,8:sum' --separater ' ' --decimal-width 0 |
                csvtk -tH cut -f 1,2,4,3 > {output.summary:q}
        """)


rule extract_flanking:
    input:
        bed=rules.summary_tnpb.output.details,
        genome="results/{genome}/genomic.fna",
        index="results/{genome}/genomic.fna.fai",
    output:
        fasta=directory("results/{genome}/flanking"),
    params:
        extends=1500,
    conda:
        "../env/bedtools_csvtk.yaml"
    shell:
        textwrap.dedent(r"""
        mkdir -p {output.fasta:q}
        bedtools slop -i {input.bed:q} -g {input.index:q} -b {params.extends} |
            bedtools getfasta -s -fi {input.genome:q} -bed - -bedOut |
            awk -v 'FS=\t' -v 'prefix={output.fasta:q}' '{{
                print ">"$4, "cluster="$9, "loc="$1":"$2"-"$3","$6, "tnpa="$8 >(prefix "/" $9 ".fna")
                print $10 >(prefix "/" $9 ".fna")
            }}'
        """)
