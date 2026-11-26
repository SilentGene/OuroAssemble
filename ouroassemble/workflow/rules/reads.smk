# 01 — read filtering: drop reads shorter than min_read_len


rule filter_reads:
    input:
        READS,
    output:
        FILTERED_READS,
    params:
        min_len=config["min_read_len"],
    threads: config["default_threads"]
    resources:
        mem_mb=cap_mem(4000),
        runtime=240,
    log:
        "logs/filter_reads.log",
    benchmark:
        "benchmarks/filter_reads.tsv"
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit seq --threads {threads} --min-len {params.min_len} {input} "
        "-o {output} 2> {log}"


rule read_stats:
    input:
        raw=READS,
        filtered=FILTERED_READS,
    output:
        "01_reads/read_stats.tsv",
    threads: config["default_threads"]
    resources:
        mem_mb=cap_mem(4000),
        runtime=120,
    log:
        "logs/read_stats.log",
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit stats --all --tabular --threads {threads} {input.raw} {input.filtered} "
        "> {output} 2> {log}"
