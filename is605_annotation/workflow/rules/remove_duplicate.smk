import textwrap


rule merge_tnpb:
    input:
        denovo_seqs = expand(
            "results/elements_list/{tnpb}/multi_copy/{tnpb}.faa",
            tnpb=glob_wildcards("results/elements_list/{tnpb}/flanking", followlinks=True).tnpb
        ),
        known_seqs = config["dedup"]["known_seqs"]
    output:
        merged = "results/dedup/merged_tnpb.faa",
    shell:
        textwrap.dedent(r"""
            xargs cat << EOF > {output.merged:q}
                {input.denovo_seqs:q} {input.known_seqs:q}
            EOF
        """)


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
        identity=config["dedup"]["identity"],
        coverage=config["dedup"]["coverage"],
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


rule filter_known_candidates:
    input:
        table="results/dedup/merged_cluster.tsv",
        known_seqs = config["dedup"]["known_seqs"],
        rep_seqs="results/dedup/merged_rep_seq.fasta",
    output:
        new_seqs = "results/dedup/candidates.faa",
    conda:
        "../env/seqkit_csvtk.yaml"
    shell:
        r"""
            seqkit seq -ni {input.known_seqs:q} |
                {{ csvtk -tH grep -f 2 -P - {input.table:q} || echo '#'; }} |
                cut -f 1 | sort -u |
                {{ csvtk -tH grep -f 1 -P - -v {input.table:q} || true; }} |
                cut -f 1 | sort -u |
                seqkit grep -f - {input.rep_seqs:q} > {output.new_seqs:q}
        """
