# 05 — junction validation: map HiFi reads to (rotated candidates + all other
#      contigs of the same assembly) and evaluate every circular junction


rule map_reads:
    input:
        ref="04_candidates/{asm}/validation_ref.fa",
        reads=FILTERED_READS,
    output:
        bam="05_validation/{asm}/reads.bam",
        bai="05_validation/{asm}/reads.bam.bai",
    threads: job("mapping")["threads"]
    resources:
        mem_mb=cap_mem(job("mapping")["mem_mb"]),
        runtime=job("mapping")["runtime"],
    log:
        "logs/map_reads/{asm}.log",
    benchmark:
        "benchmarks/map_reads/{asm}.tsv"
    conda:
        "../envs/mapping.yaml"
    shell:
        # -Y: soft-clip supplementary alignments so split reads are visible as clips
        "(minimap2 -ax map-hifi -Y -t {threads} {input.ref} {input.reads} "
        "| samtools sort -@ {threads} -m 1G -o {output.bam} - "
        "&& samtools index {output.bam}) > {log} 2>&1"


rule eval_junctions:
    input:
        bam="05_validation/{asm}/reads.bam",
        bai="05_validation/{asm}/reads.bam.bai",
        ref="04_candidates/{asm}/validation_ref.fa",
        tsv="04_candidates/{asm}/trim_rotate.tsv",
    output:
        tsv="05_validation/{asm}/junction_eval.tsv",
        sites="05_validation/{asm}/suspicious_sites.bed",
        igv="05_validation/{asm}/igv_batch.txt",
        sites_tsv="05_validation/{asm}/igv_sites.bed",
    params:
        script=SCRIPTS / "eval_junctions.py",
        j=config["junction"],
    resources:
        mem_mb=cap_mem(8000),
        runtime=240,
    log:
        "logs/eval_junctions/{asm}.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --bam {input.bam} --ref {input.ref} --candidates {input.tsv} "
        "--min-mapq {params.j[min_mapq]} --flank {params.j[flank]} "
        "--min-clip {params.j[min_clip]} --clip-window {params.j[clip_window]} "
        "--min-spanning {params.j[min_spanning]} --min-span-ratio {params.j[min_span_ratio]} "
        "--min-depth {params.j[min_depth]} "
        "--repeat-low-mapq-frac {params.j[repeat_low_mapq_frac]} "
        "--clip-frac {params.j[clip_frac]} --min-clip-reads {params.j[min_clip_reads]} "
        "--cov-window {params.j[cov_window]} --low-cov-frac {params.j[low_cov_frac]} "
        "--out-tsv {output.tsv} --out-sites {output.sites} "
        "--out-igv {output.igv} --out-igv-bed {output.sites_tsv} > {log} 2>&1"


rule igv_report:
    input:
        sites="05_validation/{asm}/igv_sites.bed",
        bam="05_validation/{asm}/reads.bam",
        bai="05_validation/{asm}/reads.bam.bai",
        ref="04_candidates/{asm}/validation_ref.fa",
    output:
        "05_validation/{asm}/igv_report.html",
    resources:
        mem_mb=cap_mem(8000),
        runtime=240,
    log:
        "logs/igv_report/{asm}.log",
    conda:
        "../envs/igv_reports.yaml"
    shell:
        "if [ -s {input.sites} ]; then "
        "  create_report {input.sites} --fasta {input.ref} --tracks {input.bam} "
        "  --flanking 2000 --output {output} > {log} 2>&1; "
        "else "
        "  echo '<html><body><p>No circular candidates.</p></body></html>' > {output}; "
        "fi"
