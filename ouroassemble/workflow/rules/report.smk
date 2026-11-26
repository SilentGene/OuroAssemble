# 09 — final tables, FASTAs and a self-contained HTML report


def report_inputs(wildcards):
    inp = {
        "read_stats": "01_reads/read_stats.tsv",
        "contigs": expand("03_contigs/{asm}/contigs.tsv", asm=ASSEMBLERS),
        "classification": "06_classify/classification.tsv",
        "clusters": "07_derep/clusters.tsv",
        "fasta": "06_classify/selected.fa",
    }
    if config["merge"]["enabled"]:
        inp["merge"] = "08_merge/selection.tsv"
        inp["mapping"] = "08_merge/mapping_stats.tsv"
    return inp


def report_opts(wildcards, input):
    if "merge" in input.keys():
        return f"--merge-selection {input.merge} --mapping-stats {input.mapping}"
    return ""


rule publish_merged:
    input:
        fasta="08_merge/merged_assembly.fa",
        tsv="08_merge/selection.tsv",
    output:
        fasta="final/merged_assembly.fa",
        tsv="final/merged_assembly.tsv",
    resources:
        mem_mb=cap_mem(1000),
        runtime=30,
    shell:
        "cp {input.fasta} {output.fasta} && cp {input.tsv} {output.tsv}"


rule report:
    input:
        unpack(report_inputs),
    output:
        html="final/report.html",
        tsv="final/circular_genomes.tsv",
        chromosomes="final/chromosomes.fa",
        plasmids="final/plasmids.fa",
        viruses="final/viruses.fa",
        other="final/other_circular.fa",
        pending="final/pending_review.fa",
    params:
        script=SCRIPTS / "make_report.py",
        sample=config["sample"],
        outdir="final",
        opts=report_opts,
        config_json=lambda w: json.dumps(
            {k: config[k] for k in ("min_read_len", "assemblers", "junction", "classify", "derep", "merge")}
        ),
    resources:
        mem_mb=cap_mem(8000),
        runtime=60,
    log:
        "logs/report.log",
    conda:
        "../envs/pyscripts.yaml"
    shell:
        "python {params.script} --sample {params.sample} --read-stats {input.read_stats} "
        "--contigs {input.contigs} --classification {input.classification} "
        "--clusters {input.clusters} --fasta {input.fasta} "
        "--config-json {params.config_json:q} {params.opts} --outdir {params.outdir} > {log} 2>&1"
