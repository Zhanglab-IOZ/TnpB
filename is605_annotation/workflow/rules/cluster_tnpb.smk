import textwrap


rule cluster_tnpb_within_genome:
    input:
        fasta="results/genomes/{genome}/tnpb.faa",
    output:
        rep="results/genomes/{genome}/tnpb_rep_seq.fasta",
        allseqs="results/genomes/{genome}/tnpb_all_seqs.fasta",
        table="results/genomes/{genome}/tnpb_cluster.tsv",
    log:
        mmseqs2="results/genomes/{genome}/tnpb_mmseqs.log",
    params:
        prefix="results/genomes/{genome}/tnpb",
        identity=config["clustering"]["identity"],
        coverage=config["clustering"]["coverage"],
    conda:
        "../env/mmseqs2.yaml"
    threads: 8
    shell:
        textwrap.dedent(r"""
            if [ -s {input.fasta:q} ]; then
                tmpd=$(mktemp -d -t --dry-run mmseq_XXXXXXXX)
                mmseqs easy-cluster {input.fasta:q} {params.prefix:q} $tmpd --min-seq-id {params.identity} -c {params.coverage} --threads {threads} > {log.mmseqs2:q}
                rm -r $tmpd
            else
                touch {output:q}
            fi
        """)


rule summary_tnpb:
    input:
        tnpa="results/genomes/{genome}/tnpa.bed",
        tnpb="results/genomes/{genome}/tnpb.bed",
        tnpb_cluster=rules.cluster_tnpb_within_genome.output.table,
    output:
        details="results/genomes/{genome}/tnpb_annot.bed",
        summary="results/genomes/{genome}/summary.tsv",
    conda:
        "../env/bedtools_csvtk.yaml"
    shell:
        textwrap.dedent(r"""
            bedtools closest -a {input.tnpb:q} -b {input.tnpa:q} -d |
                cut -f 1-7,15 |
                awk -v 'FS=\t' -v 'OFS=\t' '{{ $8 = NF>=8 ? ($8<500 && $8!=-1) : 0 ; print }}' |
                {{ csvtk -tH join --left-join - {input.tnpb_cluster:q} -f "4;2" || true; }} |
                tee {output.details:q} |
                {{ csvtk -tH summary -g 9 -f '4:count,4:collapse,8:sum' --separater ' ' --decimal-width 0 || true; }} |
                awk -v "FS=\t" -v 'OFS=\t' '{{
                    print $1, $2, $4, ($4>1)?"T":"F", $3
                }}' > {output.summary:q}
        """)


rule filter_multi_copy:
    input:
        fasta="results/genomes/{genome}/tnpb.faa",
        summary="results/genomes/{genome}/summary.tsv",
    output:
        filtered=directory("results/genomes/{genome}/multi_copy"),
    conda:
        "../env/seqkit.yaml"
    shell:
        r"""
            mkdir -p {output:q}
            cat {input.summary:q} | awk -v 'FS=\t' '$4=="T"{{print $1}}' | seqkit grep -f - {input.fasta:q} |
            awk -v 'out='{output.filtered:q} '/^>/{{name = substr($0, 2); file = out "/" name ".faa"}} {{print > file}}'
        """


rule extract_flanking:
    input:
        bed=rules.summary_tnpb.output.details,
        genome="results/genomes/{genome}/genomic.fna",
        index="results/genomes/{genome}/genomic.fna.fai",
    output:
        fasta=directory("results/genomes/{genome}/flanking"),
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


rule list_elements:
    input:
        filtered=expand("results/genomes/{genome}/multi_copy", genome = config["genomes"].keys()),
        flanking = expand("results/genomes/{genome}/flanking", genome = config["genomes"].keys()),
    output:
        outdir=directory("results/elements_list")
    params:
        genomes = list(config["genomes"].keys()),
    shell:
        textwrap.dedent(r"""
        mkdir -p {output.outdir:q}

        for genome in {params.genomes:q}; do
            for element in $(ls "results/genomes/${{genome}}/multi_copy/"); do
                ln -s ../genomes/${{genome}} {output.outdir:q}/${{element%.faa}}
            done
        done
        """)
