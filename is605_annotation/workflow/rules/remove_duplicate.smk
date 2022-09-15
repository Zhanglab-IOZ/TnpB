import textwrap


rule merge_tnpb:
    input:
        denovo_seqs = expand("results/genomes/{genome}/multi_copy.faa", genome = config["input"]["genomes"].keys()),
        known_seqs = "/dev/null"
    output:
        denovo = "results/dedup/denovo_tnpb.faa",
        merged = "results/dedup/merged_tnpb.faa",
    shell:
        r"""
            cat {input.denovo_seqs:q} >{output.denovo:q}
            cat {output.denovo:q} {input.known_seqs:q} > {output.merged:q}
        """


rule cluster_tnpb_cross_genome:
    input:
        fasta="results/dedup/merged_tnpb.faa",
    output:
        rep="results/dedup/merged_rep_seq.fasta",
        allseqs="results/dedup/merged_all_seqs.fasta",
        table="results/dedup/merged_cluster.tsv",
    log:
        mmseqs2="results/dedup/merged_mmseqs.log",
    params:
        prefix="results/dedup/merged",
        identity=0.85,
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
