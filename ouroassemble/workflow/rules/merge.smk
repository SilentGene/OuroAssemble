# 08 — merge all assemblies into one non-redundant assembly (contig selection only)
#      validated circles first, then linear contigs ranked by read support and length


rule scan_contigs:
    input:
        bam="05_validation/{asm}/reads.bam",
        bai="05_validation/{asm}/reads.bam.bai",
    output:
        "05_validation/{asm}/contig_scan.tsv",
    params:
        script=SCRIPTS / "scan_contigs.py",
        j=config["junction"],
    resources:
        mem_mb=cap_mem(8000),
        runtime=480,
    log:
        "logs/scan_contigs/{asm}.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --bam {input.bam} "
        "--min-clip {params.j[min_clip]} --clip-window {params.j[clip_window]} "
        "--clip-frac {params.j[clip_frac]} --min-clip-reads {params.j[min_clip_reads]} "
        "--cov-window {params.j[cov_window]} --low-cov-frac {params.j[low_cov_frac]} "
        "--min-depth {params.j[min_depth]} --out {output} > {log} 2>&1"


rule merge_pool:
    input:
        contigs=expand("03_contigs/{asm}/contigs.tsv", asm=ASSEMBLERS),
        fasta=expand("03_contigs/{asm}/contigs.fa", asm=ASSEMBLERS),
        trimmed=expand("04_candidates/{asm}/trimmed.fa", asm=ASSEMBLERS),
        scan=expand("05_validation/{asm}/contig_scan.tsv", asm=ASSEMBLERS),
        classification="06_classify/classification.tsv",
        clusters="07_derep/clusters.tsv",
    output:
        tsv="08_merge/pool.tsv",
        fasta="08_merge/pool.fa",
    params:
        script=SCRIPTS / "merge_assemblies.py",
        min_len=config["merge"]["min_contig_len"],
        failed="--include-failed-circular" if config["merge"]["include_failed_circular"] else "",
        priority=" ".join(config["assembler_priority"]),
    resources:
        mem_mb=cap_mem(16000),
        runtime=120,
    log:
        "logs/merge_pool.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} pool --contigs {input.contigs} --fasta {input.fasta} "
        "--trimmed {input.trimmed} --scan {input.scan} "
        "--classification {input.classification} --clusters {input.clusters} "
        "--min-contig-len {params.min_len} {params.failed} --priority {params.priority} "
        "--out-tsv {output.tsv} --out-fasta {output.fasta} > {log} 2>&1"


rule merge_align:
    input:
        "08_merge/pool.fa",
    output:
        "08_merge/pool.paf",
    threads: job("merge")["threads"]
    resources:
        mem_mb=cap_mem(job("merge")["mem_mb"]),
        runtime=job("merge")["runtime"],
    log:
        "logs/merge_align.log",
    benchmark:
        "benchmarks/merge_align.tsv"
    conda:
        "../envs/mapping.yaml"
    shell:
        # all-vs-all (-X: no self hits, each pair once), base-level alignment for identity
        "minimap2 -x asm5 -c -X --secondary=yes -N 50 -t {threads} {input} {input} "
        "> {output} 2> {log}"


rule merge_select:
    input:
        pool="08_merge/pool.tsv",
        fasta="08_merge/pool.fa",
        paf="08_merge/pool.paf",
    output:
        tsv="08_merge/selection.tsv",
        fasta="08_merge/merged_assembly.fa",
    params:
        script=SCRIPTS / "merge_assemblies.py",
        m=config["merge"],
    resources:
        mem_mb=cap_mem(16000),
        runtime=120,
    log:
        "logs/merge_select.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} select --pool {input.pool} --fasta {input.fasta} "
        "--paf {input.paf} --min-identity {params.m[min_identity]} "
        "--max-covered-frac {params.m[max_covered_frac]} --min-block {params.m[min_block]} "
        "--out-tsv {output.tsv} --out-fasta {output.fasta} > {log} 2>&1"


rule map_merged:
    input:
        ref="08_merge/merged_assembly.fa",
        reads=FILTERED_READS,
    output:
        bam="08_merge/reads.bam",
        bai="08_merge/reads.bam.bai",
    threads: job("mapping")["threads"]
    resources:
        mem_mb=cap_mem(job("mapping")["mem_mb"]),
        runtime=job("mapping")["runtime"],
    log:
        "logs/map_merged.log",
    benchmark:
        "benchmarks/map_merged.tsv"
    conda:
        "../envs/mapping.yaml"
    shell:
        "(minimap2 -ax map-hifi -Y -t {threads} {input.ref} {input.reads} "
        "| samtools sort -@ {threads} -m 1G -o {output.bam} - "
        "&& samtools index {output.bam}) > {log} 2>&1"


rule mapping_stats:
    input:
        bams=expand("05_validation/{asm}/reads.bam", asm=ASSEMBLERS) + ["08_merge/reads.bam"],
    output:
        "08_merge/mapping_stats.tsv",
    params:
        names=ASSEMBLERS + ["merged"],
    threads: config["default_threads"]
    resources:
        mem_mb=cap_mem(4000),
        runtime=240,
    log:
        "logs/mapping_stats.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "(printf 'assembly\\treads\\treads_mapped\\tpct_reads_mapped\\n' > {output}; "
        "names=({params.names}); bams=({input.bams}); "
        "for i in ${{!bams[@]}}; do "
        "  total=$(samtools view -@ {threads} -c -F 0x900 ${{bams[$i]}}); "
        "  mapped=$(samtools view -@ {threads} -c -F 0x904 ${{bams[$i]}}); "
        "  awk -v n=${{names[$i]}} -v t=$total -v m=$mapped "
        "    'BEGIN{{printf \"%s\\t%d\\t%d\\t%.2f\\n\", n, t, m, (t ? 100*m/t : 0)}}' >> {output}; "
        "done) > {log} 2>&1"
