# 07 — dereplicate circular genomes recovered by several assemblers


rule skani_triangle:
    input:
        # selected sequences >= derep.small_len (written by prepare_classification)
        "06_classify/derep_large.fa",
    output:
        "07_derep/skani_edges.tsv",
    threads: config["derep"]["threads"]
    resources:
        mem_mb=cap_mem(16000),
        runtime=240,
    log:
        "logs/skani_triangle.log",
    conda:
        "../envs/skani.yaml"
    shell:
        "n=$(grep -c '^>' {input} || true); "
        "if [ \"$n\" -ge 2 ]; then "
        "  skani triangle -i {input} -E -t {threads} -o {output} > {log} 2>&1; "
        "else "
        "  printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n' > {output}; "
        "  echo 'fewer than 2 sequences to compare' > {log}; "
        "fi"


rule dereplicate:
    input:
        classification="06_classify/classification.tsv",
        fasta="06_classify/selected.fa",
        skani="07_derep/skani_edges.tsv",
    output:
        "07_derep/clusters.tsv",
    params:
        script=SCRIPTS / "dereplicate.py",
        d=config["derep"],
        priority=" ".join(config["assembler_priority"]),
    resources:
        mem_mb=cap_mem(8000),
        runtime=120,
    log:
        "logs/dereplicate.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --classification {input.classification} "
        "--fasta {input.fasta} --skani {input.skani} "
        "--min-ani {params.d[min_ani]} --min-af {params.d[min_af]} "
        "--small-len {params.d[small_len]} --priority {params.priority} "
        "--out {output} > {log} 2>&1"
